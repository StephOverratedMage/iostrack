import os
import io
import re
import sys
import yaml
import random
import requests
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def clean_accession(raw_id):
    """Extracts clean base accession from NCBI, ENA, BV-BRC, or local FASTA headers."""
    cleaned = str(raw_id).strip().strip("|")
    if "|" in cleaned:
        cleaned = cleaned.split("|")[-1]
    match = re.search(r"([A-Z]{1,6}_?[0-9]{4,})", cleaned)
    if match:
        return match.group(1)
    # Preserve BV-BRC numeric genome/contig IDs (e.g. 2697049.18423)
    if re.match(r"^\d+\.\d+", cleaned):
        return cleaned
    return cleaned.split(".")[0]


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

BVBRC_SEQ_URL = "https://www.bv-brc.org/api/genome_sequence/"
BVBRC_HEADERS = {
    "Accept": "application/json",
    "Content-Type": "application/x-www-form-urlencoded",
}
BVBRC_SELECT = "select(accession,genome_id,description,sequence,length,release_date,date_inserted)"


def clean_meta_str(val):
    """Returns empty string for NaN/null/placeholder metadata values."""
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return ""
    text = str(val).strip()
    return "" if text.lower() in MISSING_TOKENS else text

SENTINEL_DATES = {"1969-12-31", "1970-01-01", "1900-01-01", "0000-00-00", "0001-01-01"} #Required for sanitize_date

def sanitize_date(val, is_db_timestamp=False):
    """Rejects Unix epoch-0 / placeholder dates (and pre-1980 DB insertion timestamps)."""
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
        if yr < min_yr or yr > 2030:
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


def parse_content_range_total(header_val):
    """Parses total hit count from HTTP Content-Range header (e.g., 'items 0-0/123456')."""
    if not header_val:
        return 0
    m = re.search(r"/(\d+)", str(header_val))
    return int(m.group(1)) if m else 0


def bvbrc_query_sequences(rql_query, timeout=60):
    """Posts an RQL query to BV-BRC genome_sequence and returns (status_code, items_list, total_count)."""
    try:
        resp = requests.post(BVBRC_SEQ_URL, data=rql_query, headers=BVBRC_HEADERS, timeout=timeout)
        total = parse_content_range_total(resp.headers.get("Content-Range", ""))
        if resp.status_code == 200:
            data = resp.json()
            return 200, (data if isinstance(data, list) else []), total
        return resp.status_code, [], total
    except Exception as e:
        print(f"⚠️ Warning: BV-BRC request failed ({e})")
        return 0, [], 0


def interleave_item_buckets(buckets):
    """Interleaves lists of BV-BRC item dicts while deduplicating on cleaned accession."""
    out = []
    seen = set()
    max_len = max((len(b) for b in buckets), default=0)
    for idx in range(max_len):
        for b in buckets:
            if idx < len(b):
                acc = clean_accession(b[idx].get("accession", ""))
                if acc and acc not in seen:
                    seen.add(acc)
                    out.append(b[idx])
    return out


def build_bvbrc_candidate_pool(bvbrc_query, min_len, max_len, max_ref_len, target_quota, selection_strategy, random_seed):
    """Constructs a deterministic BV-BRC candidate item list for all 7 selection strategies."""
    base_filter = f"{bvbrc_query}&ge(length,{min_len})&le(length,{max_len})"
    is_two_stage = selection_strategy in (
        "Newest + Best Alignment (Two-Stage)",
        "Fewest SNPs (Two-Stage)",
    )
    db_fetch_limit = max(target_quota * 2, 60) if selection_strategy == "Fewest SNPs (Two-Stage)" else (
        target_quota * 2 if is_two_stage else target_quota
    )
    pool_limit = max(db_fetch_limit * 5, 250)

    if selection_strategy in ("Newest First (Deterministic)", "Newest + Best Alignment (Two-Stage)"):
        q = f"{base_filter}&{BVBRC_SELECT}&sort(-date_inserted,+accession)&limit({pool_limit})"
        print(f"🌍 Querying BV-BRC API ({selection_strategy}): {q}")
        status, items, _ = bvbrc_query_sequences(q)
        if status != 200:
            status, items, _ = bvbrc_query_sequences(f"{base_filter}&{BVBRC_SELECT}&limit({pool_limit})")
        items.sort(
            key=lambda x: (
                str(x.get("date_inserted") or x.get("release_date") or ""),
                clean_accession(x.get("accession", "")),
            ),
            reverse=True,
        )
        return items, db_fetch_limit

    if selection_strategy == "Oldest First (Deterministic)":
        q = f"{base_filter}&gt(date_inserted,1980-01-01T00%3A00%3A00.000Z)&{BVBRC_SELECT}&sort(+date_inserted,+accession)&limit({pool_limit})"
        print(f"🕰️ Querying BV-BRC Server-Side Sort ({selection_strategy}): {q}")
        status, items, _ = bvbrc_query_sequences(q)
        if status != 200 or not items:
            print("⚠️ Retrying BV-BRC oldest query without server-side date filter...")
            status, items, _ = bvbrc_query_sequences(
                f"{base_filter}&{BVBRC_SELECT}&sort(+date_inserted,+accession)&limit({pool_limit})"
            )
            if status != 200 or not items:
                status, items, _ = bvbrc_query_sequences(
                    f"{base_filter}&{BVBRC_SELECT}&limit({pool_limit})"
                )
            valid_items = [
                it for it in items
                if not str(it.get("date_inserted") or "").startswith(("1969", "1970", "1900", "0000"))
            ]
            items = valid_items or items

        items.sort(
            key=lambda x: (
                str(x.get("date_inserted") or x.get("release_date") or "9999"),
                clean_accession(x.get("accession", "")),
            )
        )
        return items, db_fetch_limit

    if selection_strategy == "Longest Sequences First":
        q = f"{base_filter}&{BVBRC_SELECT}&sort(-length,+accession)&limit({pool_limit})"
        print(f"📏 Querying BV-BRC Server-Side Sort ({selection_strategy}): {q}")
        status, items, _ = bvbrc_query_sequences(q)
        if status != 200:
            status, items, _ = bvbrc_query_sequences(f"{base_filter}&{BVBRC_SELECT}&limit({pool_limit})")
        items.sort(
            key=lambda x: (
                -int(x.get("length") or len(str(x.get("sequence", "")))),
                clean_accession(x.get("accession", "")),
            )
        )
        return items, db_fetch_limit

    if selection_strategy == "Closest to Reference Length":
        print(f"📏 Querying progressive BV-BRC length bands centered around {max_ref_len:,} bp...")
        collected = []
        seen = set()
        for delta in (50, 250, 1000, None):
            if delta is None:
                lo_b, hi_b = min_len, max_len
            else:
                lo_b = max(min_len, max_ref_len - delta)
                hi_b = min(max_len, max_ref_len + delta)
            band_q = f"{bvbrc_query}&ge(length,{lo_b})&le(length,{hi_b})&{BVBRC_SELECT}&sort(+accession)&limit({pool_limit})"
            _, items, _ = bvbrc_query_sequences(band_q)
            for it in items:
                acc = clean_accession(it.get("accession", ""))
                if acc and acc not in seen:
                    seen.add(acc)
                    collected.append(it)
            if len(collected) >= pool_limit:
                break
        collected.sort(
            key=lambda x: (
                abs(int(x.get("length") or len(str(x.get("sequence", "")))) - max_ref_len),
                clean_accession(x.get("accession", "")),
            )
        )
        return collected, db_fetch_limit

    # Probe total archive count cleanly via limit(1)
    _, _, total_count = bvbrc_query_sequences(f"{base_filter}&select(accession)&sort(+accession)&limit(1)", timeout=30)

    if selection_strategy == "Fewest SNPs (Two-Stage)":
        print(f"🎯 Stage 1 ({selection_strategy}): Sampling oldest + cross-archive BV-BRC pool (total: {total_count:,})...")
        buckets = []
        _, oldest_items, _ = bvbrc_query_sequences(
            f"{base_filter}&gt(date_inserted,1980-01-01T00%3A00%3A00.000Z)&{BVBRC_SELECT}&sort(+date_inserted,+accession)&limit({max(db_fetch_limit, 80)})"
        )
        if oldest_items:
            buckets.append(oldest_items)

        safe_span = min(max(0, total_count - 40), 80000)
        if safe_span > 40:
            for frac in (0.15, 0.40, 0.70, 0.95):
                off = int(safe_span * frac)
                _, slice_items, _ = bvbrc_query_sequences(
                    f"{base_filter}&{BVBRC_SELECT}&sort(+accession)&limit(40,{off})"
                )
                if slice_items:
                    buckets.append(slice_items)
        else:
            _, acc_items, _ = bvbrc_query_sequences(
                f"{base_filter}&{BVBRC_SELECT}&sort(+accession)&limit({pool_limit})"
            )
            if acc_items:
                buckets.append(acc_items)

        return interleave_item_buckets(buckets), db_fetch_limit

    # Random Sample (Reproducible Seed): multi-offset sampling across +accession, -accession, +date_inserted, -date_inserted
    rng = random.Random(random_seed)
    print(f"🎲 Sampling BV-BRC across archive with reproducible seed={random_seed} (total: {total_count:,})...")
    if total_count <= 500:
        _, items, _ = bvbrc_query_sequences(f"{base_filter}&{BVBRC_SELECT}&sort(+accession)&limit(500)")
        items.sort(key=lambda x: clean_accession(x.get("accession", "")))
        rng.shuffle(items)
        return items, db_fetch_limit

    slice_size = max(15, (target_quota * 4) // 8)
    safe_span = max(0, min(total_count - slice_size, 100000))
    sort_modes = ["+accession", "-accession", "+date_inserted", "-date_inserted"]
    buckets = []

    for idx in range(8):
        sort_expr = sort_modes[idx % len(sort_modes)]
        off = rng.randint(0, safe_span) if safe_span > 0 else 0
        q = f"{base_filter}&{BVBRC_SELECT}&sort({sort_expr})&limit({slice_size},{off})"
        _, slice_items, _ = bvbrc_query_sequences(q)
        if slice_items:
            slice_items.sort(key=lambda x: clean_accession(x.get("accession", "")))
            rng.shuffle(slice_items)
            buckets.append(slice_items)

    combined = interleave_item_buckets(buckets)
    if not combined:
        _, combined, _ = bvbrc_query_sequences(f"{base_filter}&{BVBRC_SELECT}&sort(+accession)&limit({pool_limit})")
        combined.sort(key=lambda x: clean_accession(x.get("accession", "")))
        rng.shuffle(combined)
    return combined, db_fetch_limit


def fetch_bvbrc_genome_metadata(genome_ids):
    """Fetches structured isolation_country and collection_date from BV-BRC's genome endpoint."""
    unique_ids = sorted({str(gid).strip() for gid in genome_ids if clean_meta_str(gid)})
    if not unique_ids:
        return {}

    meta_by_gid = {}
    url = "https://www.bv-brc.org/api/genome/"

    chunk_size = 50
    for i in range(0, len(unique_ids), chunk_size):
        chunk = unique_ids[i : i + chunk_size]
        query = (
            f"in(genome_id,({','.join(chunk)}))"
            f"&select(genome_id,isolation_country,collection_date,collection_year,date_inserted,date_modified, geographic_location)"
            f"&limit({len(chunk) + 10})"
        )
        try:
            resp = requests.post(url, data=query, headers=BVBRC_HEADERS, timeout=30)
            if resp.status_code == 200:
                for row in resp.json():
                    gid = str(row.get("genome_id", "")).strip()
                    if gid:
                        meta_by_gid[gid] = row
        except Exception as e:
            print(f"⚠️ Warning: BV-BRC genome metadata lookup skipped ({e})")

    return meta_by_gid


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


def main():
    with open("config/config.yaml", "r") as f:
        config = yaml.safe_load(f)

    os.makedirs("results/raw", exist_ok=True)
    Entrez.email = config.get("ncbi_email") or "test@example.com"

    # 1. Fetch Reference Genome from NCBI (required for rich GenBank CDS features)
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

    # 2. Build BV-BRC Query
    bvbrc_query = config.get("bvbrc_query", "eq(taxon_id,11250)")
    selection_strategy = config.get("selection_strategy", "Newest First (Deterministic)")
    random_seed = int(config.get("random_seed", 42))
    target_quota = int(config.get("max_sequences", 20))
    min_len = int(config.get("min_length", 10000))
    max_len = max(50000, int(max_ref_len * 1.5))
    max_n_percent = float(config.get("max_ambiguous_n_percent", 1.0))
    search_mode = config.get("search_mode", "By Organism & Keyword")

    seen_accessions = set(original_ref_accs + [base_ref])
    local_records, local_meta = load_local_fasta_isolates(
        config, seen_accessions, min_len, max_n_percent
    )

    db_records = []
    db_meta_rows = []
    meta_cols = ["Accession", "Title", "Length_bp", "UpdateDate", "Collection_Date", "Country"]

    if search_mode != "Local FASTA Only":
        if search_mode == "By Accession List":
            full_query = f"{bvbrc_query}&{BVBRC_SELECT}"
            print(f"🌍 Querying BV-BRC API (By Accession List): {full_query}")
            status, records_data, _ = bvbrc_query_sequences(full_query)
            if status != 200 and not local_records:
                print(f"❌ BV-BRC API Error {status}", file=sys.stderr)
                sys.exit(1)
            db_fetch_limit = max(target_quota, len(records_data))
        else:
            records_data, db_fetch_limit = build_bvbrc_candidate_pool(
                bvbrc_query=bvbrc_query,
                min_len=min_len,
                max_len=max_len,
                max_ref_len=max_ref_len,
                target_quota=target_quota,
                selection_strategy=selection_strategy,
                random_seed=random_seed,
            )

        if not records_data and not local_records:
            print("❌ BV-BRC returned NO sequences for query.", file=sys.stderr)
            sys.exit(1)

        print(f"🧬 Actively hunting for up to {db_fetch_limit} pristine sequence(s) from {len(records_data)} BV-BRC candidate(s)...")

        accepted_items = []
        for item in records_data:
            raw_acc = item.get("accession", "")
            clean_id = clean_accession(raw_acc)
            seq_str = str(item.get("sequence", "")).upper()

            if not seq_str:
                continue

            n_count = seq_str.count("N")
            n_percent = (n_count / len(seq_str)) * 100.0

            if clean_id not in seen_accessions and len(seq_str) >= min_len and n_percent <= max_n_percent:
                seen_accessions.add(clean_id)
                rec = SeqRecord(
                    Seq(seq_str),
                    id=clean_id,
                    name=clean_id,
                    description="",
                )
                db_records.append(rec)
                accepted_items.append((clean_id, len(seq_str), item))

            if len(db_records) >= db_fetch_limit:
                break

        genome_ids = [item.get("genome_id", "") for _, _, item in accepted_items]
        print(f"📊 Fetching BV-BRC spatiotemporal metadata for {len(accepted_items)} sequence(s)...")
        genome_meta_map = fetch_bvbrc_genome_metadata(genome_ids)

        for clean_id, seq_len, item in accepted_items:
            raw_acc = item.get("accession", clean_id)
            gid = str(item.get("genome_id", "")).strip()
            gmeta = genome_meta_map.get(gid, {})

            title = clean_meta_str(item.get("description", "")) or raw_acc
            update_date = (
                sanitize_date(gmeta.get("date_modified", ""), is_db_timestamp=True)
                or sanitize_date(item.get("date_inserted", ""), is_db_timestamp=True)
                or sanitize_date(item.get("release_date", ""), is_db_timestamp=True)
            )

            # Replaced raw_col_date & fallback_date block:
            raw_date_field = clean_meta_str(gmeta.get("collection_date", ""))
            if raw_date_field:
                # If collection_date is a sentinel like 1969-12-31, do not fall back to its derived collection_year (1969)
                raw_col_date = sanitize_date(raw_date_field, is_db_timestamp=False)
            else:
                raw_col_date = sanitize_date(gmeta.get("collection_year", ""), is_db_timestamp=False)

            fallback_date = (
                sanitize_date(item.get("release_date", ""), is_db_timestamp=True)
                or sanitize_date(item.get("date_inserted", ""), is_db_timestamp=True)
                or update_date
            )

            country, col_date = parse_spatiotemporal_fields(
                title=title,
                raw_country=clean_meta_str(gmeta.get("isolation_country", "")) or clean_meta_str(gmeta.get("geographic_location", "")),
                raw_date=raw_col_date,
                fallback_date=fallback_date,
            )

            db_meta_rows.append({
                "Accession": clean_id,
                "Title": title,
                "Length_bp": int(item.get("length") or seq_len),
                "UpdateDate": update_date,
                "Collection_Date": col_date,
                "Country": country,
            })

    valid_records = local_records + db_records
    if not valid_records:
        print("❌ No BV-BRC or local sequences passed the QC filters.", file=sys.stderr)
        sys.exit(1)

    SeqIO.write(valid_records, "results/raw/isolates.fasta", "fasta")

    all_meta = local_meta + db_meta_rows
    final_meta_df = pd.DataFrame(all_meta, columns=meta_cols) if all_meta else pd.DataFrame(columns=meta_cols)
    final_meta_df.to_csv("results/raw/metadata.csv", index=False)
    print(f"✅ Success! Saved {len(valid_records)} pristine sequence(s) and metadata.")


if __name__ == "__main__":
    main()