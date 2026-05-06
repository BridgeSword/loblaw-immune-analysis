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


PROJECT_ROOT = Path(__file__).resolve().parent
UNKNOWN_RESPONSE_LABEL = "unknown"
TRIAL_RESPONSE_QUERY = (
    "condition == 'melanoma' and treatment == 'miraclib' "
    "and sample_type == 'PBMC' and response in ['yes', 'no']"
)


st.set_page_config(
    page_title="Loblaw Bio Immune Trial",
    page_icon="LB",
    layout="wide",
    menu_items={},
)


@st.cache_data(show_spinner=False)
def load_outputs() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    if (
        not analysis.DATABASE_PATH.exists()
        or not analysis.FREQUENCY_SUMMARY_CSV_PATH.exists()
    ):
        subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "analysis.py")],
            check=True,
            cwd=PROJECT_ROOT,
        )

    with sqlite3.connect(analysis.DATABASE_PATH) as database_connection:
        frequency_summary = analysis.read_frequency_summary(database_connection)

    response_statistics = pd.read_csv(analysis.RESPONSE_STATISTICS_CSV_PATH)
    baseline_subset_summary = json.loads(
        analysis.BASELINE_SUBSET_JSON_PATH.read_text(encoding="utf-8")
    )
    return frequency_summary, response_statistics, baseline_subset_summary


frequency_summary, response_statistics, baseline_subset_summary = load_outputs()
frequency_summary["response_filter"] = frequency_summary["response"].fillna(
    UNKNOWN_RESPONSE_LABEL
)

st.markdown(
    """
    <style>
    header[data-testid="stHeader"],
    div[data-testid="stToolbar"],
    div[data-testid="stDecoration"],
    div[data-testid="stStatusWidget"],
    div[data-testid="stHeaderActionElements"],
    button[kind="header"],
    .stDeployButton,
    #MainMenu {
        display: none !important;
        visibility: hidden !important;
    }
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

miraclib_response_frequency_summary = frequency_summary.query(TRIAL_RESPONSE_QUERY)

metric_cols = st.columns(4)
metric_cols[0].metric("Samples", f"{frequency_summary['sample'].nunique():,}")
metric_cols[1].metric("Subjects", f"{frequency_summary['subject'].nunique():,}")
metric_cols[2].metric("Projects", f"{frequency_summary['project'].nunique():,}")
metric_cols[3].metric(
    "Miraclib melanoma PBMC",
    f"{miraclib_response_frequency_summary['sample'].nunique():,}",
)

with st.sidebar:
    st.header("Filters")
    projects = st.multiselect(
        "Project",
        sorted(frequency_summary["project"].unique()),
        default=sorted(frequency_summary["project"].unique()),
    )
    conditions = st.multiselect(
        "Condition",
        sorted(frequency_summary["condition"].unique()),
        default=sorted(frequency_summary["condition"].unique()),
    )
    sample_types = st.multiselect(
        "Sample type",
        sorted(frequency_summary["sample_type"].unique()),
        default=sorted(frequency_summary["sample_type"].unique()),
    )
    treatments = st.multiselect(
        "Treatment",
        sorted(frequency_summary["treatment"].unique()),
        default=sorted(frequency_summary["treatment"].unique()),
    )
    responses = st.multiselect(
        "Response",
        sorted(frequency_summary["response_filter"].unique()),
        default=sorted(frequency_summary["response_filter"].unique()),
    )
    time_points = st.multiselect(
        "Days from treatment start",
        sorted(frequency_summary["time_from_treatment_start"].unique()),
        default=sorted(frequency_summary["time_from_treatment_start"].unique()),
    )

filtered_frequency_summary = frequency_summary[
    frequency_summary["project"].isin(projects)
    & frequency_summary["condition"].isin(conditions)
    & frequency_summary["sample_type"].isin(sample_types)
    & frequency_summary["treatment"].isin(treatments)
    & frequency_summary["response_filter"].isin(responses)
    & frequency_summary["time_from_treatment_start"].isin(time_points)
].copy()

tab_overview, tab_response, tab_subset = st.tabs(
    ["Frequency Overview", "Response Statistics", "Baseline Subset"]
)

with tab_overview:
    left, right = st.columns([1.35, 1])
    with left:
        frequency_boxplot = px.box(
            filtered_frequency_summary,
            x="population",
            y="percentage",
            color="population",
            points="outliers",
            labels={
                "population": "Population",
                "percentage": "Relative frequency (%)",
            },
        )
        frequency_boxplot.update_layout(
            showlegend=False, height=470, margin=dict(l=10, r=10, t=30, b=10)
        )
        st.plotly_chart(frequency_boxplot, use_container_width=True)
    with right:
        mean_frequency_by_population = (
            filtered_frequency_summary.groupby("population", as_index=False)["percentage"]
            .mean()
            .rename(columns={"percentage": "mean_percentage"})
            .sort_values("mean_percentage", ascending=False)
        )
        st.dataframe(
            mean_frequency_by_population,
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
        filtered_frequency_summary[
            ["sample", "total_count", "population", "count", "percentage"]
        ],
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
        miraclib_response_frequency_summary,
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
        response_statistics.sort_values("q_value_bh"),
        use_container_width=True,
        hide_index=True,
        column_config={
            "p_value": st.column_config.NumberColumn("p value", format="%.3e"),
            "q_value_bh": st.column_config.NumberColumn("BH q value", format="%.3e"),
            "significant_at_fdr_0_05": "FDR < 0.05",
        },
    )

with tab_subset:
    baseline_subset = pd.read_csv(analysis.BASELINE_SUBSET_CSV_PATH)
    c1, c2, c3 = st.columns(3)
    c1.metric("Baseline samples", f"{baseline_subset['sample'].nunique():,}")
    c2.metric(
        "Responder subjects",
        baseline_subset_summary["subjects_by_response"].get("yes", 0),
    )
    c3.metric(
        "Male responder mean B cells",
        f"{baseline_subset_summary['melanoma_male_responder_avg_b_cells_time0']:.2f}",
    )

    left, middle, right = st.columns(3)
    left.bar_chart(
        pd.Series(baseline_subset_summary["sample_count_by_project"]), height=260
    )
    middle.bar_chart(
        pd.Series(baseline_subset_summary["subjects_by_response"]), height=260
    )
    right.bar_chart(pd.Series(baseline_subset_summary["subjects_by_sex"]), height=260)

    st.dataframe(baseline_subset, use_container_width=True, hide_index=True, height=360)
