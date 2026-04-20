from __future__ import annotations

import subprocess
from pathlib import Path

from shared.run_store import (
    append_log,
    artifact_key,
    ensure_run_dirs,
    list_artifacts,
    load_metadata,
    run_prefix,
    save_metadata,
)
from shared.storage import get_storage
from worker.parser import (
    GroupRecord,
    SampleRecord,
    build_count_matrix,
    build_count_matrix_from_gene_matrix,
    parse_group_sheet,
    parse_sample_sheet,
)


def _select_records(records: list[SampleRecord], group_a: str | None, group_b: str | None) -> tuple[list[SampleRecord], str, str]:
    groups = sorted({r.group for r in records})
    if group_a is None or group_b is None:
        if len(groups) != 2:
            raise ValueError(
                f"Sample sheet must contain exactly two groups for default mode; found: {groups}"
            )
        group_a, group_b = groups[0], groups[1]

    if group_a == group_b:
        raise ValueError("group_a and group_b must be different")

    selected = [r for r in records if r.group in {group_a, group_b}]
    counts = {group_a: 0, group_b: 0}
    for rec in selected:
        counts[rec.group] += 1

    if counts[group_a] < 2 or counts[group_b] < 2:
        raise ValueError(f"Need at least two replicates per group. Observed: {counts}")

    return selected, group_a, group_b


def _select_group_records(records: list[GroupRecord], group_a: str | None, group_b: str | None) -> tuple[list[GroupRecord], str, str]:
    groups = sorted({r.group for r in records})
    if group_a is None or group_b is None:
        if len(groups) != 2:
            raise ValueError(
                f"Sample sheet must contain exactly two groups for default mode; found: {groups}"
            )
        group_a, group_b = groups[0], groups[1]

    if group_a == group_b:
        raise ValueError("group_a and group_b must be different")

    selected = [r for r in records if r.group in {group_a, group_b}]
    counts = {group_a: 0, group_b: 0}
    for rec in selected:
        counts[rec.group] += 1

    if counts[group_a] < 2 or counts[group_b] < 2:
        raise ValueError(f"Need at least two replicates per group. Observed: {counts}")

    return selected, group_a, group_b


def _download_inputs_for_run(run_id: str, uploads_dir: Path) -> None:
    st = get_storage()
    rp = run_prefix(run_id)
    files = st.list_files(f"{rp}/uploads")
    for name in files:
        st.download_to(f"{rp}/uploads/{name}", uploads_dir / name)


def _upload_artifacts_for_run(run_id: str, artifacts_dir: Path) -> list[str]:
    st = get_storage()
    names = sorted([p.name for p in artifacts_dir.glob("*") if p.is_file()])
    for name in names:
        st.upload_from(artifacts_dir / name, artifact_key(run_id, name))
    return names


def run_pipeline(run_id: str, group_a: str | None = None, group_b: str | None = None) -> None:
    meta = load_metadata(run_id)
    dirs = ensure_run_dirs(run_id)
    meta.status = "running"
    meta.error = None
    save_metadata(meta)

    try:
        input_type = meta.input_type
        append_log(run_id, f"Preparing input files for input_type={input_type}")
        _download_inputs_for_run(run_id, dirs["uploads"])

        sample_sheet = dirs["uploads"] / "sample_sheet.csv"
        append_log(run_id, f"Parsing sample sheet for input_type={input_type}")

        count_csv = dirs["work"] / "count_matrix.csv"

        if input_type == "rsem_expected_counts":
            records = parse_sample_sheet(sample_sheet)
            selected, group_a, group_b = _select_records(records, group_a, group_b)
            append_log(run_id, f"Building count matrix from RSEM files for {len(selected)} samples")
            n_samples, n_genes = build_count_matrix(selected, dirs["uploads"], count_csv)
        elif input_type in {"raw_counts_matrix", "salmon_kallisto_estimated_counts"}:
            records = parse_group_sheet(sample_sheet)
            selected, group_a, group_b = _select_group_records(records, group_a, group_b)
            matrix_path = dirs["uploads"] / "count_matrix.csv"
            append_log(run_id, f"Building count matrix from uploaded matrix for {len(selected)} samples")
            n_samples, n_genes = build_count_matrix_from_gene_matrix(selected, matrix_path, count_csv)
        elif input_type == "tpm_fpkm_only":
            raise ValueError(
                "TPM/FPKM-only input is not supported for DESeq2. Upload raw/expected/estimated counts instead."
            )
        else:
            raise ValueError(f"Unsupported input_type: {input_type}")

        meta.selected_groups = [group_a, group_b]
        save_metadata(meta)

        append_log(run_id, f"Count matrix ready: samples={n_samples}, genes={n_genes}")

        script = Path(__file__).resolve().parent / "r" / "deseq2_two_group.R"
        cmd = [
            "Rscript",
            str(script),
            "--count-matrix",
            str(count_csv),
            "--group-a",
            group_a,
            "--group-b",
            group_b,
            "--out-dir",
            str(dirs["artifacts"]),
        ]
        append_log(run_id, "Running DESeq2")
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.stdout:
            append_log(run_id, proc.stdout.strip())
        if proc.stderr:
            append_log(run_id, proc.stderr.strip())
        if proc.returncode != 0:
            raise RuntimeError(f"DESeq2 script failed with code {proc.returncode}")

        artifact_names = _upload_artifacts_for_run(run_id, dirs["artifacts"])
        meta.status = "completed"
        meta.artifacts = artifact_names
        save_metadata(meta)
        append_log(run_id, "Run completed")
    except Exception as exc:
        meta = load_metadata(run_id)
        meta.status = "failed"
        meta.error = str(exc)
        save_metadata(meta)
        append_log(run_id, f"ERROR: {exc}")
        raise
