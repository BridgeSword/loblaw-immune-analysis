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

from load_data import DATABASE_PATH as SQLITE_DATABASE_PATH

PROJECT_ROOT = Path(__file__).resolve().parent
DATABASE_PATH = SQLITE_DATABASE_PATH
OUTPUT_DIRECTORY = PROJECT_ROOT / "outputs"
PLOT_DIRECTORY = PROJECT_ROOT / "plots"
FREQUENCY_SUMMARY_CSV_PATH = OUTPUT_DIRECTORY / "frequency_summary.csv"
RESPONSE_STATISTICS_CSV_PATH = OUTPUT_DIRECTORY / "miraclib_response_statistics.csv"
BASELINE_SUBSET_CSV_PATH = OUTPUT_DIRECTORY / "baseline_melanoma_miraclib_pbmc_samples.csv"
BASELINE_SUBSET_JSON_PATH = OUTPUT_DIRECTORY / "baseline_subset_summary.json"
RESPONSE_BOXPLOT_PATH = PLOT_DIRECTORY / "miraclib_response_boxplot.png"

REQUIRED_FREQUENCY_COLUMNS = ["sample", "total_count", "population", "count", "percentage"]
RESPONSE_LABEL_ORDER = ["yes", "no"]
RESPONSE_COLOR_MAP = {"yes": "#0f766e", "no": "#b42318"}
MELANOMA_MIRACLIB_PBMC_SQL_FILTER = """
    sub.condition = 'melanoma'
    AND sub.treatment = 'miraclib'
    AND s.sample_type = 'PBMC'
"""
MELANOMA_MIRACLIB_PBMC_QUERY_FILTER = (
    "condition == 'melanoma' and treatment == 'miraclib' and sample_type == 'PBMC'"
)


def ensure_database() -> None:
    if DATABASE_PATH.exists():
        return
    subprocess.run(
        [sys.executable, str(PROJECT_ROOT / "load_data.py")],
        check=True,
        cwd=PROJECT_ROOT,
    )


def read_frequency_summary(database_connection: sqlite3.Connection) -> pd.DataFrame:
    frequency_summary_query = """
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
    return pd.read_sql_query(frequency_summary_query, database_connection)


def export_required_frequency_table(frequency_summary: pd.DataFrame) -> None:
    frequency_summary[REQUIRED_FREQUENCY_COLUMNS].to_csv(
        FREQUENCY_SUMMARY_CSV_PATH, index=False
    )


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    number_of_tests = len(p_values)
    ordered_p_values = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted_p_values = [1.0] * number_of_tests
    running_min = 1.0

    for rank, (original_index, p_value) in reversed(
        list(enumerate(ordered_p_values, start=1))
    ):
        adjusted_value = min(running_min, p_value * number_of_tests / rank)
        running_min = adjusted_value
        adjusted_p_values[original_index] = min(adjusted_value, 1.0)

    return adjusted_p_values


def filter_melanoma_miraclib_pbmc_samples(
    frequency_summary: pd.DataFrame,
) -> pd.DataFrame:
    return frequency_summary.query(MELANOMA_MIRACLIB_PBMC_QUERY_FILTER).copy()


def compare_response_groups(
    response_frequency_summary: pd.DataFrame,
) -> pd.DataFrame:
    statistical_result_rows = []
    for population_name, population_frequency_table in response_frequency_summary.groupby(
        "population", sort=True
    ):
        responder_percentages = population_frequency_table.loc[
            population_frequency_table["response"] == "yes", "percentage"
        ].astype(float)
        non_responder_percentages = population_frequency_table.loc[
            population_frequency_table["response"] == "no", "percentage"
        ].astype(float)

        if responder_percentages.empty or non_responder_percentages.empty:
            p_value = float("nan")
            statistic = float("nan")
        else:
            test_result = mannwhitneyu(
                responder_percentages,
                non_responder_percentages,
                alternative="two-sided",
            )
            p_value = float(test_result.pvalue)
            statistic = float(test_result.statistic)

        statistical_result_rows.append(
            {
                "population": population_name,
                "n_responder_samples": int(responder_percentages.shape[0]),
                "n_non_responder_samples": int(non_responder_percentages.shape[0]),
                "mean_responder_percentage": round(float(responder_percentages.mean()), 4),
                "mean_non_responder_percentage": round(
                    float(non_responder_percentages.mean()), 4
                ),
                "median_responder_percentage": round(
                    float(responder_percentages.median()), 4
                ),
                "median_non_responder_percentage": round(
                    float(non_responder_percentages.median()), 4
                ),
                "mean_difference_responder_minus_non": round(
                    float(
                        responder_percentages.mean()
                        - non_responder_percentages.mean()
                    ),
                    4,
                ),
                "mannwhitneyu_statistic": statistic,
                "p_value": p_value,
            }
        )

    response_statistics = pd.DataFrame(statistical_result_rows)
    response_statistics["significant_at_p_0_05"] = response_statistics["p_value"] < 0.05
    response_statistics["q_value_bh"] = benjamini_hochberg(
        response_statistics["p_value"].fillna(1.0).tolist()
    )
    response_statistics["significant_at_fdr_0_05"] = (
        response_statistics["q_value_bh"] < 0.05
    )
    return response_statistics


def run_response_statistics(
    frequency_summary: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    response_frequency_summary = filter_melanoma_miraclib_pbmc_samples(frequency_summary)
    response_statistics = compare_response_groups(response_frequency_summary)
    response_statistics.to_csv(RESPONSE_STATISTICS_CSV_PATH, index=False)
    return response_frequency_summary, response_statistics


def plot_response_boxplot(response_frequency_summary: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 7))
    axes = sns.boxplot(
        data=response_frequency_summary,
        x="population",
        y="percentage",
        hue="response",
        hue_order=RESPONSE_LABEL_ORDER,
        palette=RESPONSE_COLOR_MAP,
        showfliers=False,
    )
    sns.stripplot(
        data=response_frequency_summary,
        x="population",
        y="percentage",
        hue="response",
        hue_order=RESPONSE_LABEL_ORDER,
        dodge=True,
        alpha=0.18,
        size=2,
        palette=RESPONSE_COLOR_MAP,
        ax=axes,
    )
    legend_handles, _ = axes.get_legend_handles_labels()
    axes.legend(legend_handles[:2], ["Responder", "Non-responder"], title="Response")
    axes.set_title("Miraclib response comparison in melanoma PBMC samples")
    axes.set_xlabel("Immune cell population")
    axes.set_ylabel("Relative frequency (%)")
    plt.tight_layout()
    plt.savefig(RESPONSE_BOXPLOT_PATH, dpi=180)
    plt.close()


def run_baseline_subset_analysis(
    database_connection: sqlite3.Connection,
) -> dict[str, object]:
    baseline_subset_query = f"""
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
        {MELANOMA_MIRACLIB_PBMC_SQL_FILTER}
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
    baseline_subset = pd.read_sql_query(baseline_subset_query, database_connection)
    baseline_subset.to_csv(BASELINE_SUBSET_CSV_PATH, index=False)

    male_responder_average_b_cell_count = baseline_subset.query(
        "sex == 'M' and response == 'yes'"
    )["b_cell"].mean()
    baseline_subset_summary = {
        "filter": {
            "condition": "melanoma",
            "sample_type": "PBMC",
            "treatment": "miraclib",
            "time_from_treatment_start": 0,
            "note": "quintazide is mentioned for completeness, but it is not a treatment value in this dataset.",
        },
        "sample_count_by_project": baseline_subset.groupby("project")["sample"]
        .nunique()
        .astype(int)
        .to_dict(),
        "subjects_by_response": baseline_subset.groupby("response")["subject"]
        .nunique()
        .astype(int)
        .to_dict(),
        "subjects_by_sex": baseline_subset.groupby("sex")["subject"]
        .nunique()
        .astype(int)
        .to_dict(),
        "melanoma_male_responder_avg_b_cells_time0": round(
            float(male_responder_average_b_cell_count), 2
        ),
    }

    BASELINE_SUBSET_JSON_PATH.write_text(
        json.dumps(baseline_subset_summary, indent=2), encoding="utf-8"
    )
    return baseline_subset_summary


def main() -> None:
    OUTPUT_DIRECTORY.mkdir(exist_ok=True)
    PLOT_DIRECTORY.mkdir(exist_ok=True)
    ensure_database()

    with sqlite3.connect(DATABASE_PATH) as database_connection:
        frequency_summary = read_frequency_summary(database_connection)
        export_required_frequency_table(frequency_summary)
        response_frequency_summary, response_statistics = run_response_statistics(
            frequency_summary
        )
        plot_response_boxplot(response_frequency_summary)
        baseline_subset_summary = run_baseline_subset_analysis(database_connection)

    fdr_significant_populations = response_statistics.loc[
        response_statistics["significant_at_fdr_0_05"], "population"
    ].tolist()
    print(f"Wrote frequency table: {FREQUENCY_SUMMARY_CSV_PATH}")
    print(f"Wrote response statistics: {RESPONSE_STATISTICS_CSV_PATH}")
    print(f"Wrote response boxplot: {RESPONSE_BOXPLOT_PATH}")
    print(f"Wrote baseline subset table: {BASELINE_SUBSET_CSV_PATH}")
    print(f"Wrote baseline subset summary: {BASELINE_SUBSET_JSON_PATH}")
    print(f"Significant populations at FDR 0.05: {fdr_significant_populations or 'none'}")
    print(
        "Melanoma male responder average B cells at time 0: "
        f"{baseline_subset_summary['melanoma_male_responder_avg_b_cells_time0']:.2f}"
    )


if __name__ == "__main__":
    main()
