"""Interactive Streamlit dashboard for the Loblaw Bio immune trial analysis."""

from __future__ import annotations

import json
import sqlite3
import subprocess
import sys
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

import analysis


ROOT = Path(__file__).resolve().parent


st.set_page_config(
    page_title="Loblaw Bio Immune Trial",
    page_icon="LB",
    layout="wide",
)


@st.cache_data(show_spinner=False)
def load_outputs() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    if not analysis.DB_PATH.exists() or not analysis.SUMMARY_CSV.exists():
        subprocess.run([sys.executable, str(ROOT / "analysis.py")], check=True, cwd=ROOT)

    with sqlite3.connect(analysis.DB_PATH) as conn:
        summary = analysis.read_frequency_summary(conn)

    stats = pd.read_csv(analysis.STATS_CSV)
    subset_summary = json.loads(analysis.SUBSET_JSON.read_text(encoding="utf-8"))
    return summary, stats, subset_summary


summary_df, stats_df, subset_summary = load_outputs()
summary_df["response_filter"] = summary_df["response"].fillna("unknown")

st.markdown(
    """
    <style>
    .block-container {padding-top: 1.4rem; padding-bottom: 2rem;}
    div[data-testid="stMetric"] {
        border: 1px solid #d8dde3;
        border-radius: 8px;
        padding: 0.85rem 1rem;
        background: #ffffff;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("Loblaw Bio Immune Trial Dashboard")

filtered_trial = summary_df.query(
    "condition == 'melanoma' and treatment == 'miraclib' and sample_type == 'PBMC' and response in ['yes', 'no']"
)

metric_cols = st.columns(4)
metric_cols[0].metric("Samples", f"{summary_df['sample'].nunique():,}")
metric_cols[1].metric("Subjects", f"{summary_df['subject'].nunique():,}")
metric_cols[2].metric("Projects", f"{summary_df['project'].nunique():,}")
metric_cols[3].metric("Miraclib melanoma PBMC", f"{filtered_trial['sample'].nunique():,}")

with st.sidebar:
    st.header("Filters")
    projects = st.multiselect(
        "Project",
        sorted(summary_df["project"].unique()),
        default=sorted(summary_df["project"].unique()),
    )
    conditions = st.multiselect(
        "Condition",
        sorted(summary_df["condition"].unique()),
        default=sorted(summary_df["condition"].unique()),
    )
    sample_types = st.multiselect(
        "Sample type",
        sorted(summary_df["sample_type"].unique()),
        default=sorted(summary_df["sample_type"].unique()),
    )
    treatments = st.multiselect(
        "Treatment",
        sorted(summary_df["treatment"].unique()),
        default=sorted(summary_df["treatment"].unique()),
    )
    responses = st.multiselect(
        "Response",
        sorted(summary_df["response_filter"].unique()),
        default=sorted(summary_df["response_filter"].unique()),
    )
    time_points = st.multiselect(
        "Days from treatment start",
        sorted(summary_df["time_from_treatment_start"].unique()),
        default=sorted(summary_df["time_from_treatment_start"].unique()),
    )

filtered = summary_df[
    summary_df["project"].isin(projects)
    & summary_df["condition"].isin(conditions)
    & summary_df["sample_type"].isin(sample_types)
    & summary_df["treatment"].isin(treatments)
    & summary_df["response_filter"].isin(responses)
    & summary_df["time_from_treatment_start"].isin(time_points)
].copy()

tab_overview, tab_response, tab_subset = st.tabs(
    ["Frequency Overview", "Response Statistics", "Baseline Subset"]
)

with tab_overview:
    left, right = st.columns([1.35, 1])
    with left:
        fig = px.box(
            filtered,
            x="population",
            y="percentage",
            color="population",
            points="outliers",
            labels={
                "population": "Population",
                "percentage": "Relative frequency (%)",
            },
        )
        fig.update_layout(showlegend=False, height=470, margin=dict(l=10, r=10, t=30, b=10))
        st.plotly_chart(fig, use_container_width=True)
    with right:
        totals = (
            filtered.groupby("population", as_index=False)["percentage"]
            .mean()
            .rename(columns={"percentage": "mean_percentage"})
            .sort_values("mean_percentage", ascending=False)
        )
        st.dataframe(
            totals,
            use_container_width=True,
            hide_index=True,
            column_config={
                "population": "Population",
                "mean_percentage": st.column_config.NumberColumn(
                    "Mean frequency (%)", format="%.2f"
                ),
            },
        )

    st.dataframe(
        filtered[["sample", "total_count", "population", "count", "percentage"]],
        use_container_width=True,
        hide_index=True,
        height=330,
        column_config={
            "percentage": st.column_config.NumberColumn("Percentage", format="%.2f")
        },
    )

with tab_response:
    st.subheader("Melanoma PBMC samples treated with miraclib")
    response_fig = px.box(
        filtered_trial,
        x="population",
        y="percentage",
        color="response",
        points="all",
        category_orders={"response": ["yes", "no"]},
        color_discrete_map={"yes": "#0f766e", "no": "#b42318"},
        labels={
            "population": "Population",
            "percentage": "Relative frequency (%)",
            "response": "Response",
        },
    )
    response_fig.update_layout(height=520, margin=dict(l=10, r=10, t=25, b=10))
    st.plotly_chart(response_fig, use_container_width=True)

    st.dataframe(
        stats_df.sort_values("q_value_bh"),
        use_container_width=True,
        hide_index=True,
        column_config={
            "p_value": st.column_config.NumberColumn("p value", format="%.3e"),
            "q_value_bh": st.column_config.NumberColumn("BH q value", format="%.3e"),
            "significant_at_fdr_0_05": "FDR < 0.05",
        },
    )

with tab_subset:
    subset = pd.read_csv(analysis.SUBSET_CSV)
    c1, c2, c3 = st.columns(3)
    c1.metric("Baseline samples", f"{subset['sample'].nunique():,}")
    c2.metric("Responder subjects", subset_summary["subjects_by_response"].get("yes", 0))
    c3.metric(
        "Male responder mean B cells",
        f"{subset_summary['melanoma_male_responder_avg_b_cells_time0']:.2f}",
    )

    left, middle, right = st.columns(3)
    left.bar_chart(pd.Series(subset_summary["sample_count_by_project"]), height=260)
    middle.bar_chart(pd.Series(subset_summary["subjects_by_response"]), height=260)
    right.bar_chart(pd.Series(subset_summary["subjects_by_sex"]), height=260)

    st.dataframe(subset, use_container_width=True, hide_index=True, height=360)
