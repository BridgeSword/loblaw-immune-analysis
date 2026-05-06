"""Run Loblaw Bio immune cell analyses and export reproducible outputs."""

from __future__ import annotations

import json
import sqlite3
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "immune_trial.db"
OUTPUT_DIR = ROOT / "outputs"
PLOT_DIR = ROOT / "plots"
SUMMARY_CSV = OUTPUT_DIR / "frequency_summary.csv"
STATS_CSV = OUTPUT_DIR / "miraclib_response_statistics.csv"
SUBSET_CSV = OUTPUT_DIR / "baseline_melanoma_miraclib_pbmc_samples.csv"
SUBSET_JSON = OUTPUT_DIR / "baseline_subset_summary.json"
BOXPLOT = PLOT_DIR / "miraclib_response_boxplot.png"

FILTER_CLAUSE = """
    sub.condition = 'melanoma'
    AND sub.treatment = 'miraclib'
    AND s.sample_type = 'PBMC'
"""


def ensure_database() -> None:
    if DB_PATH.exists():
        return
    subprocess.run([sys.executable, str(ROOT / "load_data.py")], check=True, cwd=ROOT)


def read_frequency_summary(conn: sqlite3.Connection) -> pd.DataFrame:
    query = """
    WITH sample_totals AS (
        SELECT sample_id, SUM(cell_count) AS total_count
        FROM cell_counts
        GROUP BY sample_id
    )
    SELECT
        s.sample_code AS sample,
        sample_totals.total_count AS total_count,
        p.population_name AS population,
        cc.cell_count AS count,
        ROUND(100.0 * cc.cell_count / sample_totals.total_count, 4) AS percentage,
        pr.project_name AS project,
        sub.subject_code AS subject,
        sub.condition,
        sub.age,
        sub.sex,
        sub.treatment,
        sub.response,
        s.sample_type,
        s.time_from_treatment_start
    FROM cell_counts AS cc
    JOIN populations AS p ON p.population_id = cc.population_id
    JOIN samples AS s ON s.sample_id = cc.sample_id
    JOIN sample_totals ON sample_totals.sample_id = s.sample_id
    JOIN subjects AS sub ON sub.subject_id = s.subject_id
    JOIN projects AS pr ON pr.project_id = sub.project_id
    ORDER BY s.sample_code, p.population_name
    """
    return pd.read_sql_query(query, conn)


def export_required_frequency_table(summary: pd.DataFrame) -> None:
    required_columns = ["sample", "total_count", "population", "count", "percentage"]
    summary[required_columns].to_csv(SUMMARY_CSV, index=False)


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    n = len(p_values)
    ordered = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [1.0] * n
    running_min = 1.0

    for rank, (original_index, p_value) in reversed(list(enumerate(ordered, start=1))):
        candidate = min(running_min, p_value * n / rank)
        running_min = candidate
        adjusted[original_index] = min(candidate, 1.0)

    return adjusted


def run_response_statistics(summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    filtered = summary.query(
        "condition == 'melanoma' and treatment == 'miraclib' and sample_type == 'PBMC'"
    ).copy()

    rows = []
    for population, population_df in filtered.groupby("population", sort=True):
        responders = population_df.loc[
            population_df["response"] == "yes", "percentage"
        ].astype(float)
        non_responders = population_df.loc[
            population_df["response"] == "no", "percentage"
        ].astype(float)

        if responders.empty or non_responders.empty:
            p_value = float("nan")
            statistic = float("nan")
        else:
            test = mannwhitneyu(responders, non_responders, alternative="two-sided")
            p_value = float(test.pvalue)
            statistic = float(test.statistic)

        rows.append(
            {
                "population": population,
                "n_responder_samples": int(responders.shape[0]),
                "n_non_responder_samples": int(non_responders.shape[0]),
                "mean_responder_percentage": round(float(responders.mean()), 4),
                "mean_non_responder_percentage": round(float(non_responders.mean()), 4),
                "median_responder_percentage": round(float(responders.median()), 4),
                "median_non_responder_percentage": round(
                    float(non_responders.median()), 4
                ),
                "mean_difference_responder_minus_non": round(
                    float(responders.mean() - non_responders.mean()), 4
                ),
                "mannwhitneyu_statistic": statistic,
                "p_value": p_value,
            }
        )

    stats = pd.DataFrame(rows)
    stats["q_value_bh"] = benjamini_hochberg(stats["p_value"].fillna(1.0).tolist())
    stats["significant_at_fdr_0_05"] = stats["q_value_bh"] < 0.05
    stats.to_csv(STATS_CSV, index=False)
    return filtered, stats


def plot_response_boxplot(filtered: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 7))
    ax = sns.boxplot(
        data=filtered,
        x="population",
        y="percentage",
        hue="response",
        hue_order=["yes", "no"],
        palette={"yes": "#0f766e", "no": "#b42318"},
        showfliers=False,
    )
    sns.stripplot(
        data=filtered,
        x="population",
        y="percentage",
        hue="response",
        hue_order=["yes", "no"],
        dodge=True,
        alpha=0.18,
        size=2,
        palette={"yes": "#0f766e", "no": "#b42318"},
        ax=ax,
    )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], ["Responder", "Non-responder"], title="Response")
    ax.set_title("Miraclib response comparison in melanoma PBMC samples")
    ax.set_xlabel("Immune cell population")
    ax.set_ylabel("Relative frequency (%)")
    plt.tight_layout()
    plt.savefig(BOXPLOT, dpi=180)
    plt.close()


def run_subset_analysis(conn: sqlite3.Connection) -> dict[str, object]:
    query = f"""
    SELECT
        pr.project_name AS project,
        sub.subject_code AS subject,
        sub.condition,
        sub.sex,
        sub.treatment,
        sub.response,
        s.sample_code AS sample,
        s.sample_type,
        s.time_from_treatment_start,
        MAX(CASE WHEN p.population_name = 'b_cell' THEN cc.cell_count END) AS b_cell
    FROM samples AS s
    JOIN subjects AS sub ON sub.subject_id = s.subject_id
    JOIN projects AS pr ON pr.project_id = sub.project_id
    JOIN cell_counts AS cc ON cc.sample_id = s.sample_id
    JOIN populations AS p ON p.population_id = cc.population_id
    WHERE
        {FILTER_CLAUSE}
        AND s.time_from_treatment_start = 0
    GROUP BY
        pr.project_name,
        sub.subject_code,
        sub.condition,
        sub.sex,
        sub.treatment,
        sub.response,
        s.sample_code,
        s.sample_type,
        s.time_from_treatment_start
    ORDER BY pr.project_name, s.sample_code
    """
    subset = pd.read_sql_query(query, conn)
    subset.to_csv(SUBSET_CSV, index=False)

    responders_male_avg = subset.query("sex == 'M' and response == 'yes'")[
        "b_cell"
    ].mean()
    summary = {
        "filter": {
            "condition": "melanoma",
            "sample_type": "PBMC",
            "treatment": "miraclib",
            "time_from_treatment_start": 0,
            "note": "quintazide is mentioned for completeness, but it is not a treatment value in this dataset.",
        },
        "sample_count_by_project": subset.groupby("project")["sample"]
        .nunique()
        .astype(int)
        .to_dict(),
        "subjects_by_response": subset.groupby("response")["subject"]
        .nunique()
        .astype(int)
        .to_dict(),
        "subjects_by_sex": subset.groupby("sex")["subject"]
        .nunique()
        .astype(int)
        .to_dict(),
        "melanoma_male_responder_avg_b_cells_time0": round(
            float(responders_male_avg), 2
        ),
    }

    SUBSET_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    PLOT_DIR.mkdir(exist_ok=True)
    ensure_database()

    with sqlite3.connect(DB_PATH) as conn:
        summary = read_frequency_summary(conn)
        export_required_frequency_table(summary)
        filtered, stats = run_response_statistics(summary)
        plot_response_boxplot(filtered)
        subset_summary = run_subset_analysis(conn)

    significant = stats.loc[stats["significant_at_fdr_0_05"], "population"].tolist()
    print(f"Wrote frequency table: {SUMMARY_CSV}")
    print(f"Wrote response statistics: {STATS_CSV}")
    print(f"Wrote response boxplot: {BOXPLOT}")
    print(f"Wrote baseline subset table: {SUBSET_CSV}")
    print(f"Wrote baseline subset summary: {SUBSET_JSON}")
    print(f"Significant populations at FDR 0.05: {significant or 'none'}")
    print(
        "Melanoma male responder average B cells at time 0: "
        f"{subset_summary['melanoma_male_responder_avg_b_cells_time0']:.2f}"
    )


if __name__ == "__main__":
    main()
