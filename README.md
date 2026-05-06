# Loblaw Bio Immune Cell Analysis

This repository analyzes immune cell population counts from `data/cell-count.csv` for Bob Loblaw's clinical trial. It builds a SQLite database, exports frequency and statistical analysis outputs, and provides an interactive Streamlit dashboard.

Dashboard link after starting the server: [http://localhost:8501](http://localhost:8501)

## Quick Start

```bash
make setup
make pipeline
make dashboard
```

The grader can run the same three targets in GitHub Codespaces. `python load_data.py` also works directly from the repository root and creates `immune_trial.db` in the root directory.

## Outputs

Running `make pipeline` creates:

- `immune_trial.db`: SQLite database in the repository root.
- `outputs/frequency_summary.csv`: required Part 2 table with `sample`, `total_count`, `population`, `count`, and `percentage`.
- `outputs/miraclib_response_statistics.csv`: Mann-Whitney U response comparison for melanoma PBMC samples treated with miraclib.
- `plots/miraclib_response_boxplot.png`: response vs non-response boxplot for every immune cell population.
- `outputs/baseline_melanoma_miraclib_pbmc_samples.csv`: Part 4 baseline sample subset.
- `outputs/baseline_subset_summary.json`: Part 4 project, response, sex, and average B-cell summary.

The treatment name `quintazide` is mentioned here for completeness because it appears in the assignment prompt, but it is not present as a treatment value in this dataset; the requested analyses use `miraclib`.

## Database Schema

The SQLite schema is normalized into five tables:

- `projects`: one row per project, keyed by `project_id`.
- `subjects`: patient-level metadata, including project, condition, age, sex, treatment, and response.
- `samples`: sample-level metadata, including sample code, sample type, subject, and time from treatment start.
- `populations`: one row per immune population name.
- `cell_counts`: long-format fact table with one count per sample and population.

This design avoids storing repeated metadata on every cell-count row while preserving the natural relationships in the trial: projects have subjects, subjects have samples, and samples have measured cell populations. It also scales cleanly if there are hundreds of projects, thousands of samples, more time points, or new immune cell populations. New populations can be inserted into `populations` without changing the count table shape, and analytics can filter by indexed subject or sample metadata before aggregating the fact table.

## Analysis Design

`load_data.py` is intentionally in the repository root because the submission requires direct execution with `python load_data.py`. It recreates the database each run, validates the CSV columns, loads reference tables, and inserts counts in normalized long format.

`analysis.py` owns the reproducible pipeline. It reads from SQLite, generates the Part 2 frequency table, performs the Part 3 response comparison, makes the static boxplot, and exports the Part 4 baseline subset summaries.

`dashboard.py` is a Streamlit app that reads the same database and generated outputs. It provides interactive filters for the full frequency table, an interactive response comparison, and the baseline subset view Bob can use for early treatment exploration.

## Statistical Method

For Part 3, the program filters to:

- `condition = melanoma`
- `treatment = miraclib`
- `sample_type = PBMC`

For each immune population, it compares relative frequency percentages between responders (`response = yes`) and non-responders (`response = no`) using a two-sided Mann-Whitney U test. The output also includes Benjamini-Hochberg adjusted q-values and flags populations significant at FDR 0.05.

## Part 4 Query Result

The baseline subset is melanoma PBMC samples at `time_from_treatment_start = 0` from subjects treated with miraclib. The summary JSON includes:

- sample counts by project
- responder and non-responder subject counts
- male and female subject counts
- average B-cell count for melanoma male responders at time 0, rounded to two decimals
