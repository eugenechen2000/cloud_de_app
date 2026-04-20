from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

REQUIRED_SAMPLE_SHEET_COLUMNS = {"sample_id", "file_name", "group"}
REQUIRED_GROUP_SHEET_COLUMNS = {"sample_id", "group"}


@dataclass
class SampleRecord:
    sample_id: str
    file_name: str
    group: str


@dataclass
class GroupRecord:
    sample_id: str
    group: str


def parse_sample_sheet(path: Path) -> list[SampleRecord]:
    df = pd.read_csv(path)
    missing = REQUIRED_SAMPLE_SHEET_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"Sample sheet missing columns: {sorted(missing)}")

    rows: list[SampleRecord] = []
    for _, row in df.iterrows():
        rows.append(
            SampleRecord(
                sample_id=str(row["sample_id"]).strip(),
                file_name=str(row["file_name"]).strip(),
                group=str(row["group"]).strip(),
            )
        )

    if not rows:
        raise ValueError("Sample sheet has no rows")

    return rows


def parse_group_sheet(path: Path) -> list[GroupRecord]:
    df = pd.read_csv(path)
    missing = REQUIRED_GROUP_SHEET_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"Sample sheet missing columns: {sorted(missing)}")

    rows: list[GroupRecord] = []
    for _, row in df.iterrows():
        rows.append(
            GroupRecord(
                sample_id=str(row["sample_id"]).strip(),
                group=str(row["group"]).strip(),
            )
        )
    if not rows:
        raise ValueError("Sample sheet has no rows")
    return rows


def read_expected_count(path: Path) -> pd.Series:
    if not path.exists():
        raise FileNotFoundError(f"Missing RSEM file: {path}")
    df = pd.read_csv(path, sep="	", usecols=["gene_id", "expected_count"])
    genes = df[df["gene_id"].astype(str).str.startswith("ENSG")].copy()
    genes["gene_id"] = genes["gene_id"].astype(str)
    return genes.drop_duplicates("gene_id").set_index("gene_id")["expected_count"].astype(float)


def build_count_matrix(records: list[SampleRecord], uploads_dir: Path, output_csv: Path) -> tuple[int, int]:
    first_path = uploads_dir / records[0].file_name
    first = read_expected_count(first_path)
    gene_order = first.index.tolist()

    matrix_rows: list[dict[str, float | str]] = []
    for rec in records:
        series = read_expected_count(uploads_dir / rec.file_name)
        aligned = series.reindex(gene_order).fillna(0.0)
        row: dict[str, float | str] = {
            "sample_id": rec.sample_id,
            "file_name": rec.file_name,
            "group": rec.group,
        }
        row.update({gene: float(val) for gene, val in aligned.items()})
        matrix_rows.append(row)

    out_df = pd.DataFrame(matrix_rows)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_csv, index=False)
    return len(matrix_rows), len(gene_order)


def build_count_matrix_from_gene_matrix(
    records: list[GroupRecord],
    count_matrix_path: Path,
    output_csv: Path,
) -> tuple[int, int]:
    if not count_matrix_path.exists():
        raise FileNotFoundError(f"Missing count matrix file: {count_matrix_path}")

    df = pd.read_csv(count_matrix_path)
    if df.shape[1] < 3:
        raise ValueError("Count matrix must have gene_id column + >=2 sample columns")

    gene_col = str(df.columns[0])
    sample_cols = [str(c) for c in df.columns[1:]]

    wanted = [r.sample_id for r in records]
    missing = sorted(set(wanted) - set(sample_cols))
    if missing:
        raise ValueError(f"Sample IDs missing from count matrix columns: {missing}")

    genes = df[gene_col].astype(str).tolist()

    matrix_rows: list[dict[str, float | str]] = []
    group_map = {r.sample_id: r.group for r in records}
    for sid in wanted:
        vals = pd.to_numeric(df[sid], errors="coerce").fillna(0.0).tolist()
        row: dict[str, float | str] = {
            "sample_id": sid,
            "file_name": sid,
            "group": group_map[sid],
        }
        row.update({g: float(v) for g, v in zip(genes, vals)})
        matrix_rows.append(row)

    out_df = pd.DataFrame(matrix_rows)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_csv, index=False)
    return len(matrix_rows), len(genes)
