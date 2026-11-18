import os
import sys
import yaml
import subprocess
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
import time
import re

MACRO_CONFIG_PATH = "config/macro_config.yaml"
MACRO_RESULTS_PATH = "results/macro/macro_metadata.csv"

MACRO_STRATEGIES = [
    "Newest First (Deterministic)",
    "Oldest First (Deterministic)",
    "Longest Sequences First",
    "Random Sample (Reproducible Seed)",
]

MACRO_WIDGET_STATE_KEYS = (
    "macro_filter_countries",
    "macro_filter_dates",
    "macro_drill_box",
    "macro_row_selection",
)


def reset_macro_view_state():
    """Clears cached Streamlit filter/export widget state so fresh runs never show stale values."""
    for key in MACRO_WIDGET_STATE_KEYS:
        st.session_state.pop(key, None)


def flush_macro_results():
    """Deletes saved macro metadata CSV and resets cached macro UI state."""
    removed = False
    if os.path.exists(MACRO_RESULTS_PATH):
        os.remove(MACRO_RESULTS_PATH)
        removed = True
    reset_macro_view_state()
    st.session_state.pop("macro_last_logs", None)
    return removed

def parse_row_selection(spec, max_len):
    """
    Parses a human-friendly row selection string (e.g., '0-5, 10, 156-255', 'all')
    into (selected_indices, warnings_list) for a 0-indexed table of length max_len.
    Ranges are inclusive (e.g., 0-5 selects rows 0, 1, 2, 3, 4, 5).
    """
    if max_len <= 0:
        return [], []

    raw = str(spec or "").strip()
    if not raw:
        return [], []

    if raw.lower() == "all":
        return list(range(max_len)), []

    indices = []
    seen = set()
    warnings = []

    tokens = [t.strip() for t in re.split(r"[,;]+", raw) if t.strip()]
    for tok in tokens:
        m_range = re.fullmatch(r"(\d+)\s*[-:]\s*(\d+)", tok)
        if m_range:
            start, end = int(m_range.group(1)), int(m_range.group(2))
            if start > end:
                start, end = end, start
            if start >= max_len:
                warnings.append(f"Range '{tok}' is out of bounds (valid rows: 0–{max_len - 1}).")
                continue
            if end >= max_len:
                warnings.append(f"Range '{tok}' clamped to max row {max_len - 1}.")
                end = max_len - 1
            for idx in range(start, end + 1):
                if idx not in seen:
                    seen.add(idx)
                    indices.append(idx)
        elif re.fullmatch(r"\d+", tok):
            idx = int(tok)
            if idx >= max_len:
                warnings.append(f"Row '{tok}' is out of bounds (valid rows: 0–{max_len - 1}).")
                continue
            if idx not in seen:
                seen.add(idx)
                indices.append(idx)
        else:
            warnings.append(f"Ignored invalid token '{tok}'. Use numbers or ranges like '0-5, 10'.")

    return indices, warnings


def load_macro_config():
    defaults = {
        "data_source": "NCBI",
        "organism_query": "Severe acute respiratory syndrome coronavirus 2",
        "max_records": 1000,
        "min_length": 10000,
        "max_length": 500000,
        "selection_strategy": "Random Sample (Reproducible Seed)",
        "random_seed": 42,
        "ncbi_email": "test@example.com",
    }
    if os.path.exists(MACRO_CONFIG_PATH):
        try:
            with open(MACRO_CONFIG_PATH, "r") as f:
                loaded = yaml.safe_load(f) or {}
            defaults.update(loaded)
        except Exception:
            pass
    elif os.path.exists("config/config.yaml"):
        try:
            with open("config/config.yaml", "r") as f:
                deep_cfg = yaml.safe_load(f) or {}
            if deep_cfg.get("ncbi_email"):
                defaults["ncbi_email"] = deep_cfg["ncbi_email"]
        except Exception:
            pass
    return defaults


def build_choropleth_figure(filtered_df, use_log_scale=True):
    """Builds a Plotly world choropleth map from Country + ISO3 columns."""
    mapped = filtered_df[
        filtered_df["ISO3"].fillna("").astype(str).str.strip().ne("")
        & filtered_df["Country"].fillna("Unknown").astype(str).ne("Unknown")
    ]
    if mapped.empty:
        return None

    geo_counts = (
        mapped.groupby(["ISO3", "Country"], as_index=False)
        .size()
        .rename(columns={"size": "Isolate_Count"})
        .sort_values("Isolate_Count", ascending=False)
    )

    if use_log_scale:
        geo_counts["Log10_Count"] = np.log10(geo_counts["Isolate_Count"].clip(lower=1))
        max_pow = max(1, int(np.ceil(geo_counts["Log10_Count"].max())))
        tick_vals = list(range(0, max_pow + 1))
        tick_text = [f"{10 ** p:,}" for p in tick_vals]

        fig = px.choropleth(
            geo_counts,
            locations="ISO3",
            locationmode="ISO-3",
            color="Log10_Count",
            hover_name="Country",
            hover_data={"Isolate_Count": ":,", "ISO3": True, "Log10_Count": False},
            color_continuous_scale="YlOrRd",
            title="Global Surveillance Density (Log10 Scale)",
        )
        fig.update_layout(
            coloraxis_colorbar=dict(
                title="Isolates",
                tickvals=tick_vals,
                ticktext=tick_text,
            )
        )
    else:
        fig = px.choropleth(
            geo_counts,
            locations="ISO3",
            locationmode="ISO-3",
            color="Isolate_Count",
            hover_name="Country",
            hover_data={"Isolate_Count": ":,", "ISO3": True},
            color_continuous_scale="YlOrRd",
            title="Global Surveillance Density (Linear Scale)",
            labels={"Isolate_Count": "Isolates"},
        )

    fig.update_layout(
        geo=dict(
            showframe=False,
            showcoastlines=True,
            coastlinecolor="rgba(120,120,120,0.4)",
            showcountries=True,
            countrycolor="rgba(120,120,120,0.25)",
            projection_type="natural earth",
        ),
        margin=dict(l=0, r=0, t=42, b=0),
        height=480,
    )
    return fig


def run_macro_mode():
    """Renders the Macro-Epidemiological Surveillance sidebar and dashboard."""
    mcfg = load_macro_config()

    st.sidebar.markdown("---")
    st.sidebar.subheader("🌍 Macro Surveillance Settings")

    sources = ["NCBI", "ENA", "BV-BRC"]
    src_idx = sources.index(mcfg["data_source"]) if mcfg["data_source"] in sources else 0
    data_source = st.sidebar.selectbox("Database Source", sources, index=src_idx, key="macro_db")

    organism_query = st.sidebar.text_input(
        "Organism Name / TaxID / Query",
        value=str(mcfg.get("organism_query", "Severe acute respiratory syndrome coronavirus 2")),
        help="Enter an organism name (e.g., 'Severe acute respiratory syndrome coronavirus 2'), a numeric NCBI TaxID (e.g., 2697049), or a native API query.",
        key="macro_org",
    )

    max_records = st.sidebar.slider(
        "Target Metadata Records",
        min_value=100,
        max_value=5000,
        value=min(5000, max(100, int(mcfg.get("max_records", 1000)))),
        step=100,
        help="Metadata-only mode skips FASTA downloads and alignment, fetching thousands of records in seconds.",
        key="macro_n",
    )

    col_len1, col_len2 = st.sidebar.columns(2)
    with col_len1:
        min_length = st.number_input(
            "Min Length (bp)",
            min_value=100,
            max_value=10000000,
            value=int(mcfg.get("min_length", 10000)),
            step=1000,
            key="macro_min_len",
        )
    with col_len2:
        max_length = st.number_input(
            "Max Length (bp)",
            min_value=1000,
            max_value=50000000,
            value=int(mcfg.get("max_length", 500000)),
            step=10000,
            key="macro_max_len",
        )

    saved_strat = mcfg.get("selection_strategy", MACRO_STRATEGIES[-1])
    strat_idx = MACRO_STRATEGIES.index(saved_strat) if saved_strat in MACRO_STRATEGIES else 0
    selection_strategy = st.sidebar.selectbox(
        "🎯 Selection Strategy",
        options=MACRO_STRATEGIES,
        index=strat_idx,
        help="Two-stage alignment strategies require sequence downloads and are only available in Deep Genomic Mode.",
        key="macro_strategy",
    )

    random_seed = int(mcfg.get("random_seed", 42))
    if selection_strategy == "Random Sample (Reproducible Seed)":
        random_seed = st.sidebar.number_input(
            "🎲 Random Seed",
            min_value=0,
            max_value=999999,
            value=random_seed,
            step=1,
            key="macro_seed",
        )

    ncbi_email = st.sidebar.text_input(
        "NCBI Email (for Entrez API)",
        value=str(mcfg.get("ncbi_email", "test@example.com")),
        key="macro_email",
    )

    st.sidebar.markdown("---")
    auto_flush_macro = st.sidebar.checkbox(
        "🗑️ Auto-flush previous results on Run",
        value=False,
        key="macro_auto_flush",
        help="Deletes results/macro/macro_metadata.csv and resets filters right before each new fetch.",
    )
    if auto_flush_macro:
        st.sidebar.warning(
            "⚠️ **Auto-flush enabled:** Previous macro results (`results/macro/macro_metadata.csv`) "
            "will be permanently deleted every time you click Run. Export any CSV you want to keep first."
        )

    run_btn = st.sidebar.button("🌍 Run Macro Surveillance", type="primary", use_container_width=True)
    flush_btn = st.sidebar.button(
        "🧹 Flush Macro Results",
        use_container_width=True,
        help="Deletes results/macro/macro_metadata.csv and clears active filters without touching Deep Genomic results.",
    )

    st.title("🌍 Iostrack — Macro-Epidemiological Surveillance")
    st.caption(
        "High-volume, metadata-only spatiotemporal analysis (no FASTA downloads or alignment). "
        "Public sequence counts reflect genomic surveillance capacity and upload schedules alongside true incidence."
    )

    if flush_btn:
        removed = flush_macro_results()
        if removed:
            st.toast("Macro results flushed cleanly.", icon="🧹")
        else:
            st.toast("No saved macro results found to flush.", icon="ℹ️")
        st.rerun()

    if run_btn:
        os.makedirs("config", exist_ok=True)
        new_cfg = {
            "data_source": str(data_source),
            "organism_query": str(organism_query).strip(),
            "max_records": int(max_records),
            "min_length": int(min_length),
            "max_length": int(max(max_length, min_length)),
            "selection_strategy": str(selection_strategy),
            "random_seed": int(random_seed),
            "ncbi_email": str(ncbi_email).strip() or "test@example.com",
        }
        with open(MACRO_CONFIG_PATH, "w") as f:
            yaml.safe_dump(new_cfg, f, sort_keys=False)

        if auto_flush_macro:
            flush_macro_results()
        else:
            reset_macro_view_state()

        t_start = time.monotonic()
        log_lines = []
        last_pct = 0

        with st.status(
            f"🌍 Fetching up to {max_records:,} metadata records from {data_source}...",
            expanded=True,
        ) as status_box:
            progress_bar = st.progress(0, text="0% — Initializing Macro Surveillance...")
            live_stage = st.empty()

            proc = subprocess.Popen(
                [sys.executable, "-u", "workflow/scripts/fetch_macro.py", MACRO_CONFIG_PATH],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            if proc.stdout is not None:
                for raw_line in iter(proc.stdout.readline, ""):
                    line = raw_line.rstrip("\n")
                    elapsed = int(time.monotonic() - t_start)
                    if line.startswith("PROGRESS|"):
                        parts = line.split("|", 2)
                        if len(parts) == 3:
                            try:
                                pct = max(last_pct, min(100, int(parts[1])))
                            except ValueError:
                                pct = last_pct
                            msg = parts[2].strip()
                            last_pct = pct
                            progress_bar.progress(pct, text=f"{pct}% — {msg} ({elapsed}s)")
                            live_stage.caption(f"⏳ {msg}")
                            log_lines.append(f"[{pct:3d}% | {elapsed:3d}s] {msg}")
                    elif line.strip():
                        log_lines.append(line)
                        live_stage.caption(line.strip())

            proc.wait()
            total_elapsed = time.monotonic() - t_start
            st.session_state["macro_last_logs"] = "\n".join(log_lines) or "(No output produced)"

            if proc.returncode != 0:
                progress_bar.progress(last_pct, text=f"Failed at {last_pct}% ({total_elapsed:.1f}s)")
                status_box.update(
                    label=f"❌ Macro metadata fetch failed (exit code {proc.returncode})",
                    state="error",
                    expanded=True,
                )
                st.code(st.session_state["macro_last_logs"])
                return
            else:
                progress_bar.progress(100, text=f"100% — Completed in {total_elapsed:.1f}s")
                live_stage.empty()
                status_box.update(
                    label=f"✅ Macro fetch completed ({data_source} | {total_elapsed:.1f}s)",
                    state="complete",
                    expanded=False,
                )

    if st.session_state.get("macro_last_logs"):
        with st.expander("View Last Fetch Logs", expanded=False):
            st.code(st.session_state["macro_last_logs"])

    if not os.path.exists(MACRO_RESULTS_PATH):
        st.info("Configure your organism and target record count in the sidebar, then click **🌍 Run Macro Surveillance**.")
        return

    df = pd.read_csv(MACRO_RESULTS_PATH)
    if df.empty:
        st.warning("Macro metadata file is empty. Please run a new query.")
        return

    for col, default_val in [("Country", "Unknown"), ("ISO3", ""), ("Collection_Date", "Unknown"), ("Data_Source", data_source)]:
        if col not in df.columns:
            df[col] = default_val

    df["Country"] = df["Country"].fillna("Unknown").astype(str).str.strip()
    df["ISO3"] = df["ISO3"].fillna("").astype(str).str.strip()
    df["Clean_Date"] = pd.to_datetime(df["Collection_Date"], format="mixed", errors="coerce")

    # Interactive Cohort Filters
    with st.expander("🔎 Filter Macro Cohort (Country & Date Range)", expanded=False):
        f_col1, f_col2 = st.columns(2)
        all_countries = sorted(df["Country"].unique().tolist())
        with f_col1:
            selected_countries = st.multiselect(
                "Countries (leave blank for all)",
                options=all_countries,
                default=[],
                key="macro_filter_countries",
            )
        valid_dates_all = df["Clean_Date"].dropna()
        with f_col2:
            if not valid_dates_all.empty:
                min_dt = valid_dates_all.min().date()
                max_dt = valid_dates_all.max().date()
                if "macro_filter_dates" in st.session_state:
                    saved_dr = st.session_state["macro_filter_dates"]
                    if not (
                        isinstance(saved_dr, tuple)
                        and len(saved_dr) == 2
                        and min_dt <= saved_dr[0] <= max_dt
                        and min_dt <= saved_dr[1] <= max_dt
                    ):
                        st.session_state.pop("macro_filter_dates", None)
                date_range = st.date_input(
                    "Collection Date Window",
                    value=(min_dt, max_dt),
                    min_value=min_dt,
                    max_value=max_dt,
                    key="macro_filter_dates",
                )
            else:
                date_range = None

    filtered = df.copy()
    if selected_countries:
        filtered = filtered[filtered["Country"].isin(selected_countries)]
    if isinstance(date_range, tuple) and len(date_range) == 2:
        start_d, end_d = pd.Timestamp(date_range[0]), pd.Timestamp(date_range[1])
        filtered = filtered[
            filtered["Clean_Date"].isna()
            | ((filtered["Clean_Date"] >= start_d) & (filtered["Clean_Date"] <= end_d))
        ]

    filtered = filtered.reset_index(drop=True) # Reset the index for easy specifying which rows to export in the drill-down section

    # Top KPI Strip
    total_n = len(filtered)
    known_geo = filtered[filtered["Country"] != "Unknown"]
    n_countries = known_geo["Country"].nunique()
    geo_pct = (len(known_geo) / total_n * 100.0) if total_n > 0 else 0.0
    valid_dates = filtered.dropna(subset=["Clean_Date"]).sort_values("Clean_Date")
    date_span = (
        f"{valid_dates['Clean_Date'].min().strftime('%Y-%m-%d')} → {valid_dates['Clean_Date'].max().strftime('%Y-%m-%d')}"
        if not valid_dates.empty
        else "N/A"
    )

    k1, k2, k3, k4 = st.columns(4)
    k1.metric("Total Isolates", f"{total_n:,}")
    k2.metric("Countries Represented", f"{n_countries:,}")
    k3.metric("Geo-Tagged Records", f"{geo_pct:.1f}%")
    k4.metric("Date Span", date_span)

    st.markdown("---")

    # 1. Global Choropleth Map
    map_hdr1, map_hdr2 = st.columns([3, 1])
    with map_hdr1:
        st.subheader("🗺️ Global Choropleth Distribution")
    with map_hdr2:
        color_scale_mode = st.radio(
            "Color Scale",
            ["Log10 Scale", "Linear Scale"],
            horizontal=True,
            key="macro_map_scale",
        )

    fig_choro = build_choropleth_figure(filtered, use_log_scale=(color_scale_mode == "Log10 Scale"))
    if fig_choro is not None:
        st.plotly_chart(fig_choro, width="stretch")
    else:
        st.info("No mappable ISO-3 country records found in the current selection.")

    st.markdown("---")

    # 2. Standard Geotemporal Graphs (Same 3-Scale Timeline + Donut Chart as Tab 3)
    st.subheader("📈 Spatiotemporal & Geographic Breakdown")
    geo_col1, geo_col2 = st.columns(2)

    with geo_col1:
        st.markdown("##### 📈 Epidemic Timeline (Collection Dates)")
        if not valid_dates.empty:
            time_view = st.radio(
                "Timeline Scale",
                [
                    "Continuous Daily (Show Empty Days)",
                    "Active Days Only (Compact)",
                    "Monthly Bins",
                ],
                horizontal=True,
                key="macro_time_scale_radio",
            )

            if time_view == "Monthly Bins":
                timeline_counts = (
                    valid_dates.groupby(valid_dates["Clean_Date"].dt.to_period("M"))
                    .size()
                    .reset_index(name="Isolate_Count")
                )
                timeline_counts["Clean_Date"] = timeline_counts["Clean_Date"].astype(str)
                x_col, x_title, axis_type = "Clean_Date", "Collection Month", "category"

            elif time_view == "Active Days Only (Compact)":
                timeline_counts = (
                    valid_dates["Clean_Date"]
                    .dt.strftime("%Y-%m-%d")
                    .value_counts()
                    .reset_index(name="Isolate_Count")
                )
                timeline_counts.columns = ["Clean_Date", "Isolate_Count"]
                timeline_counts = timeline_counts.sort_values("Clean_Date")
                x_col, x_title, axis_type = "Clean_Date", "Active Collection Dates", "category"

            else:
                daily_series = valid_dates["Clean_Date"].dt.normalize().value_counts().sort_index()
                min_d, max_d = daily_series.index.min(), daily_series.index.max()
                if (max_d - min_d).days <= 3650:
                    full_range = pd.date_range(min_d, max_d, freq="D")
                    daily_series = daily_series.reindex(full_range, fill_value=0)
                timeline_counts = daily_series.reset_index()
                timeline_counts.columns = ["Clean_Date", "Isolate_Count"]
                x_col, x_title, axis_type = "Clean_Date", "Calendar Date (Empty Days Preserved)", "date"

            fig_time = px.bar(
                timeline_counts,
                x=x_col,
                y="Isolate_Count",
                title="Isolates Collected / Reported Over Time",
                labels={x_col: x_title, "Isolate_Count": "Number of Isolates"},
            )
            fig_time.update_xaxes(type=axis_type)
            st.plotly_chart(fig_time, width="stretch")
        else:
            st.info("No parseable collection dates available for the current cohort.")

    with geo_col2:
        st.markdown("##### 🥧 Geographic Distribution (Country)")
        country_counts = (
            filtered["Country"]
            .value_counts()
            .reset_index(name="Count")
        )
        country_counts.columns = ["Country", "Count"]

        if len(country_counts) > 12:
            pie_view = st.radio(
                "Pie Slices",
                ["Top 12 + Other", "Show All Countries"],
                horizontal=True,
                key="macro_pie_slices",
            )
            if pie_view == "Top 12 + Other":
                top12 = country_counts.head(12).copy()
                other_sum = int(country_counts.iloc[12:]["Count"].sum())
                country_plot_df = pd.concat(
                    [top12, pd.DataFrame([{"Country": f"Other ({len(country_counts) - 12} countries)", "Count": other_sum}])],
                    ignore_index=True,
                )
            else:
                country_plot_df = country_counts
        else:
            country_plot_df = country_counts

        fig_pie = px.pie(
            country_plot_df,
            names="Country",
            values="Count",
            title="Isolates by Country / Region",
            hole=0.4,
        )
        st.plotly_chart(fig_pie, width="stretch")

    st.markdown("---")

    # 3. Metadata Table & Drill-Down into Deep Genomic Mode
    filtered = filtered.reset_index(drop=True)
    exp_col1, exp_col2 = st.columns([2, 1])
    with exp_col1:
        st.markdown("##### 📋 Filtered Metadata Table")
        display_cols = [
            c
            for c in ["Accession", "Country", "Collection_Date", "UpdateDate", "Length_bp", "Title", "Data_Source"]
            if c in filtered.columns
        ]
        st.dataframe(filtered[display_cols], use_container_width=True, height=340)

    with exp_col2:
        st.markdown("##### 🧬 Drill Down to Deep Genomic Mode")
        st.caption(
            "Select row numbers from the left table to build a copyable list for "
            "**Deep Genomic Mode → By Accession List**."
        )
        max_idx = max(0, len(filtered) - 1)
        default_spec = f"0-{min(29, max_idx)}" if len(filtered) > 0 else ""

        row_spec = st.text_input(
            f"Row Selection (0–{max_idx}, or 'all')",
            value=default_spec,
            placeholder="e.g., 0-5, 10, 156-255 or all",
            help=(
                "Use inclusive ranges and single row numbers matching the index column on the left table "
                "(e.g., '0-5, 10, 156-255'). Duplicates are removed automatically."
            ),
            key="macro_row_selection",
        )

        selected_indices, parse_warnings = parse_row_selection(row_spec, len(filtered))
        for w in parse_warnings:
            st.warning(w)

        if selected_indices:
            drill_accs = (
                filtered.iloc[selected_indices]["Accession"]
                .dropna()
                .astype(str)
                .drop_duplicates()
                .tolist()
            )
        else:
            drill_accs = []

        st.caption(f"Selected **{len(drill_accs):,}** unique accession(s).")
        if len(drill_accs) > 50:
            st.caption("⚠️ *Tip: Deep Genomic alignment is fastest with ≤ 50 isolates at a time.*")

        st.text_area(
            "Copyable Accession List",
            value=", ".join(drill_accs),
            height=130,
        )
        st.download_button(
            "⬇️ Download Full Filtered CSV",
            data=filtered.drop(columns=["Clean_Date"], errors="ignore").to_csv(index=False).encode("utf-8"),
            file_name="iostrack_macro_metadata.csv",
            mime="text/csv",
            use_container_width=True,
        )