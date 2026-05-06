"""Initialize and load the Loblaw Bio immune cell SQLite database."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "data" / "cell-count.csv"
DB_PATH = ROOT / "immune_trial.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
METADATA_COLUMNS = [
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
]


SCHEMA = """
PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS populations;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id INTEGER PRIMARY KEY,
    project_name TEXT NOT NULL UNIQUE
);

CREATE TABLE subjects (
    subject_id INTEGER PRIMARY KEY,
    project_id INTEGER NOT NULL,
    subject_code TEXT NOT NULL,
    condition TEXT NOT NULL,
    age INTEGER,
    sex TEXT NOT NULL CHECK (sex IN ('M', 'F')),
    treatment TEXT NOT NULL,
    response TEXT CHECK (response IS NULL OR response IN ('yes', 'no')),
    FOREIGN KEY (project_id) REFERENCES projects(project_id),
    UNIQUE (project_id, subject_code)
);

CREATE TABLE samples (
    sample_id INTEGER PRIMARY KEY,
    sample_code TEXT NOT NULL UNIQUE,
    subject_id INTEGER NOT NULL,
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL,
    FOREIGN KEY (subject_id) REFERENCES subjects(subject_id)
);

CREATE TABLE populations (
    population_id INTEGER PRIMARY KEY,
    population_name TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id INTEGER NOT NULL,
    population_id INTEGER NOT NULL,
    cell_count INTEGER NOT NULL CHECK (cell_count >= 0),
    PRIMARY KEY (sample_id, population_id),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    FOREIGN KEY (population_id) REFERENCES populations(population_id)
);

CREATE INDEX idx_subjects_trial_filter
    ON subjects(condition, treatment, response, sex);

CREATE INDEX idx_samples_time_type
    ON samples(sample_type, time_from_treatment_start);
"""


def read_source_data() -> pd.DataFrame:
    """Read and validate the source CSV."""
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"Could not find input CSV at {CSV_PATH}")

    df = pd.read_csv(CSV_PATH)
    expected = set(METADATA_COLUMNS + POPULATIONS)
    missing = sorted(expected - set(df.columns))
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    return df[METADATA_COLUMNS + POPULATIONS].copy()


def initialize_database(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA)
    conn.executemany(
        "INSERT INTO populations (population_name) VALUES (?)",
        [(population,) for population in POPULATIONS],
    )
    conn.commit()


def load_projects(conn: sqlite3.Connection, df: pd.DataFrame) -> dict[str, int]:
    projects = sorted(df["project"].dropna().unique())
    conn.executemany(
        "INSERT INTO projects (project_name) VALUES (?)",
        [(project,) for project in projects],
    )
    conn.commit()
    return dict(conn.execute("SELECT project_name, project_id FROM projects").fetchall())


def load_subjects(
    conn: sqlite3.Connection, df: pd.DataFrame, project_ids: dict[str, int]
) -> dict[tuple[int, str], int]:
    subject_rows = (
        df[
            [
                "project",
                "subject",
                "condition",
                "age",
                "sex",
                "treatment",
                "response",
            ]
        ]
        .drop_duplicates(["project", "subject"])
        .sort_values(["project", "subject"])
    )
    records = [
        (
            project_ids[row.project],
            row.subject,
            row.condition,
            int(row.age),
            row.sex,
            row.treatment,
            None if pd.isna(row.response) else row.response,
        )
        for row in subject_rows.itertuples(index=False)
    ]
    conn.executemany(
        """
        INSERT INTO subjects
            (project_id, subject_code, condition, age, sex, treatment, response)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        records,
    )
    conn.commit()

    return {
        (project_id, subject_code): subject_id
        for subject_id, project_id, subject_code in conn.execute(
            "SELECT subject_id, project_id, subject_code FROM subjects"
        )
    }


def load_samples_and_counts(
    conn: sqlite3.Connection,
    df: pd.DataFrame,
    project_ids: dict[str, int],
    subject_ids: dict[tuple[int, str], int],
) -> None:
    population_ids = dict(
        conn.execute("SELECT population_name, population_id FROM populations").fetchall()
    )

    sample_records = []
    count_records = []
    next_sample_id = 1

    for row in df.sort_values("sample").itertuples(index=False):
        project_id = project_ids[row.project]
        subject_id = subject_ids[(project_id, row.subject)]
        sample_records.append(
            (
                next_sample_id,
                row.sample,
                subject_id,
                row.sample_type,
                int(row.time_from_treatment_start),
            )
        )
        for population in POPULATIONS:
            count_records.append(
                (
                    next_sample_id,
                    population_ids[population],
                    int(getattr(row, population)),
                )
            )
        next_sample_id += 1

    conn.executemany(
        """
        INSERT INTO samples
            (sample_id, sample_code, subject_id, sample_type, time_from_treatment_start)
        VALUES (?, ?, ?, ?, ?)
        """,
        sample_records,
    )
    conn.executemany(
        """
        INSERT INTO cell_counts (sample_id, population_id, cell_count)
        VALUES (?, ?, ?)
        """,
        count_records,
    )
    conn.commit()


def main() -> None:
    df = read_source_data()
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute("PRAGMA foreign_keys = ON")
        initialize_database(conn)
        project_ids = load_projects(conn, df)
        subject_ids = load_subjects(conn, df, project_ids)
        load_samples_and_counts(conn, df, project_ids, subject_ids)

    print(f"Loaded {len(df):,} samples into {DB_PATH.name}")


if __name__ == "__main__":
    main()
