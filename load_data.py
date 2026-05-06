"""Initialize and load the Loblaw Bio immune cell SQLite database."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
SOURCE_CSV_PATH = PROJECT_ROOT / "data" / "cell-count.csv"
DATABASE_PATH = PROJECT_ROOT / "immune_trial.db"

CELL_POPULATION_COLUMNS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
SOURCE_METADATA_COLUMNS = [
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


DATABASE_SCHEMA = """
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
    if not SOURCE_CSV_PATH.exists():
        raise FileNotFoundError(f"Could not find input CSV at {SOURCE_CSV_PATH}")

    source_data = pd.read_csv(SOURCE_CSV_PATH)
    required_columns = set(SOURCE_METADATA_COLUMNS + CELL_POPULATION_COLUMNS)
    missing_columns = sorted(required_columns - set(source_data.columns))
    if missing_columns:
        raise ValueError(f"Input CSV is missing required columns: {missing_columns}")

    return source_data[SOURCE_METADATA_COLUMNS + CELL_POPULATION_COLUMNS].copy()


def initialize_database(database_connection: sqlite3.Connection) -> None:
    database_connection.executescript(DATABASE_SCHEMA)
    database_connection.executemany(
        "INSERT INTO populations (population_name) VALUES (?)",
        [(population_name,) for population_name in CELL_POPULATION_COLUMNS],
    )
    database_connection.commit()


def load_projects(
    database_connection: sqlite3.Connection, source_data: pd.DataFrame
) -> dict[str, int]:
    project_names = sorted(source_data["project"].dropna().unique())
    database_connection.executemany(
        "INSERT INTO projects (project_name) VALUES (?)",
        [(project_name,) for project_name in project_names],
    )
    database_connection.commit()
    return dict(
        database_connection.execute(
            "SELECT project_name, project_id FROM projects"
        ).fetchall()
    )


def load_subjects(
    database_connection: sqlite3.Connection,
    source_data: pd.DataFrame,
    project_ids_by_name: dict[str, int],
) -> dict[tuple[int, str], int]:
    unique_subject_rows = (
        source_data[
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
    subject_records = [
        (
            project_ids_by_name[source_row.project],
            source_row.subject,
            source_row.condition,
            int(source_row.age),
            source_row.sex,
            source_row.treatment,
            None if pd.isna(source_row.response) else source_row.response,
        )
        for source_row in unique_subject_rows.itertuples(index=False)
    ]
    database_connection.executemany(
        """
        INSERT INTO subjects
            (project_id, subject_code, condition, age, sex, treatment, response)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        subject_records,
    )
    database_connection.commit()

    return {
        (project_id, subject_code): subject_id
        for subject_id, project_id, subject_code in database_connection.execute(
            "SELECT subject_id, project_id, subject_code FROM subjects"
        )
    }


def load_samples_and_counts(
    database_connection: sqlite3.Connection,
    source_data: pd.DataFrame,
    project_ids_by_name: dict[str, int],
    subject_ids_by_project_and_code: dict[tuple[int, str], int],
) -> None:
    population_ids_by_name = dict(
        database_connection.execute(
            "SELECT population_name, population_id FROM populations"
        ).fetchall()
    )

    sample_records = []
    cell_count_records = []
    next_sample_id = 1

    for source_row in source_data.sort_values("sample").itertuples(index=False):
        project_id = project_ids_by_name[source_row.project]
        subject_id = subject_ids_by_project_and_code[(project_id, source_row.subject)]
        sample_records.append(
            (
                next_sample_id,
                source_row.sample,
                subject_id,
                source_row.sample_type,
                int(source_row.time_from_treatment_start),
            )
        )
        for population_name in CELL_POPULATION_COLUMNS:
            cell_count_records.append(
                (
                    next_sample_id,
                    population_ids_by_name[population_name],
                    int(getattr(source_row, population_name)),
                )
            )
        next_sample_id += 1

    database_connection.executemany(
        """
        INSERT INTO samples
            (sample_id, sample_code, subject_id, sample_type, time_from_treatment_start)
        VALUES (?, ?, ?, ?, ?)
        """,
        sample_records,
    )
    database_connection.executemany(
        """
        INSERT INTO cell_counts (sample_id, population_id, cell_count)
        VALUES (?, ?, ?)
        """,
        cell_count_records,
    )
    database_connection.commit()


def main() -> None:
    source_data = read_source_data()
    with sqlite3.connect(DATABASE_PATH) as database_connection:
        database_connection.execute("PRAGMA foreign_keys = ON")
        initialize_database(database_connection)
        project_ids_by_name = load_projects(database_connection, source_data)
        subject_ids_by_project_and_code = load_subjects(
            database_connection, source_data, project_ids_by_name
        )
        load_samples_and_counts(
            database_connection,
            source_data,
            project_ids_by_name,
            subject_ids_by_project_and_code,
        )

    print(f"Loaded {len(source_data):,} samples into {DATABASE_PATH.name}")


if __name__ == "__main__":
    main()
