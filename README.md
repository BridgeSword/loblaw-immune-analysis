# Loblaw Bio Immune Cell Analysis

This repository analyzes immune cell population counts from `data/cell-count.csv` for Bob Loblaw's clinical trial. It builds a SQLite database, exports frequency and statistical analysis outputs, and provides an interactive Streamlit dashboard.

Dashboard link after starting the server: [http://localhost:8501](http://localhost:8501)

The dashboard is included because the assignment asks for an interactive dashboard to display Bob's analysis results. It is implemented as a local Streamlit app rather than a separately hosted website, so it is launched with `make dashboard` in Codespaces or locally.

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
- `outputs/miraclib_response_statistics.csv`: Mann-Whitney U response comparison for melanoma PBMC samples treated with miraclib, including raw p-values and Benjamini-Hochberg q-values.
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

## Results and Conclusions

The source file contains `10,500` biological samples from `3,500` subjects across `3` projects. Each sample has counts for five immune populations: `b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, and `monocyte`. The frequency summary expands the dataset into long format with one row per sample and immune population, so the required Part 2 output contains `52,500` rows.

For the miraclib response analysis, the program filters to melanoma PBMC samples only. Among samples with recorded response labels, this comparison includes `993` responder sample-population observations per immune population and `975` non-responder sample-population observations per immune population. The Mann-Whitney U tests show the strongest unadjusted signal for `cd4_t_cell`, where responders have a higher mean relative frequency than non-responders (`30.5378%` vs `29.9023%`, raw p-value `0.01334`). However, after Benjamini-Hochberg multiple-testing correction across the five immune populations, the `cd4_t_cell` q-value is `0.06670`, which is above the FDR 0.05 threshold.

The main statistical conclusion is therefore conservative: there is not enough evidence at FDR 0.05 to claim that any measured immune cell population has a significant relative-frequency difference between miraclib responders and non-responders in melanoma PBMC samples. The raw `cd4_t_cell` signal is still worth reporting to Bob and Yah D'yada as a hypothesis-generating trend, but it should not be presented as a statistically significant biomarker without additional validation.

For baseline early-treatment exploration, the program identifies `656` melanoma PBMC baseline samples (`time_from_treatment_start = 0`) from miraclib-treated subjects. These samples come from `prj1` and `prj3`: `384` samples from `prj1` and `272` samples from `prj3`. The baseline subset contains `331` responder subjects and `325` non-responder subjects, with `344` male subjects and `312` female subjects. Among melanoma male responders at baseline, the average B-cell count is `10401.28`.

## Statistical Method

For Part 3, the program filters to:

- `condition = melanoma`
- `treatment = miraclib`
- `sample_type = PBMC`

For each immune population, it compares relative frequency percentages between responders (`response = yes`) and non-responders (`response = no`) using a two-sided Mann-Whitney U test. The output includes both unadjusted p-value flags and Benjamini-Hochberg adjusted q-values. The FDR-adjusted result is used for the final significance conclusion because five immune populations are tested.

## Part 4 Query Result

The baseline subset is melanoma PBMC samples at `time_from_treatment_start = 0` from subjects treated with miraclib. The summary JSON includes:

- sample counts by project
- responder and non-responder subject counts
- male and female subject counts
- average B-cell count for melanoma male responders at time 0, rounded to two decimals

For the provided CSV, the melanoma male responder average B-cell count at time 0 is `10401.28`.
