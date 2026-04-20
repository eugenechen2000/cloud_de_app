from __future__ import annotations

import csv
import uuid
from datetime import datetime, timezone
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles

from backend.app.queue import get_queue
from backend.app.schemas import (
    ArtifactItem,
    ArtifactListResponse,
    CreateRunResponse,
    StartRunRequest,
    StartRunResponse,
)
from shared.config import RUN_INLINE
from shared.models import InputType, RunMetadata
from shared.run_store import (
    artifact_key,
    ensure_run_dirs,
    list_artifacts,
    load_metadata,
    run_dir,
    run_prefix,
    save_metadata,
)
from shared.storage import get_storage
from worker.jobs import process_run_job

app = FastAPI(title="Cloud Upload-to-DESeq2")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

static_dir = Path(__file__).resolve().parents[1] / "static"
app.mount("/static", StaticFiles(directory=static_dir), name="static")

COUNT_BASED_INPUT_TYPES = {
    "rsem_expected_counts",
    "raw_counts_matrix",
    "salmon_kallisto_estimated_counts",
}


def _now() -> datetime:
    return datetime.now(tz=timezone.utc)


def _parse_groups_from_bytes(sample_sheet_bytes: bytes, required_cols: set[str]) -> list[str]:
    text = sample_sheet_bytes.decode("utf-8")
    reader = csv.DictReader(text.splitlines())
    if not required_cols.issubset(set(reader.fieldnames or [])):
        raise HTTPException(
            status_code=400,
            detail=f"Sample sheet requires columns: {','.join(sorted(required_cols))}",
        )
    groups = sorted({(row.get("group") or "").strip() for row in reader if (row.get("group") or "").strip()})
    return groups


@app.get("/")
def index() -> FileResponse:
    return FileResponse(static_dir / "index.html")


@app.post("/projects/{project_id}/runs", response_model=CreateRunResponse)
def create_run(project_id: str) -> CreateRunResponse:
    run_id = str(uuid.uuid4())
    ensure_run_dirs(run_id)
    meta = RunMetadata(
        run_id=run_id,
        project_id=project_id,
        status="created",
        created_at=_now(),
        updated_at=_now(),
    )
    save_metadata(meta)
    return CreateRunResponse(
        run_id=run_id,
        project_id=project_id,
        status=meta.status,
        created_at=meta.created_at,
    )


@app.post("/runs/{run_id}/upload-rsem")
async def upload_rsem_inputs(
    run_id: str,
    sample_sheet: UploadFile = File(...),
    rsem_files: list[UploadFile] = File(...),
) -> dict:
    try:
        meta = load_metadata(run_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    st = get_storage()
    rp = run_prefix(run_id)

    sample_bytes = await sample_sheet.read()
    if not sample_bytes:
        raise HTTPException(status_code=400, detail="Sample sheet is empty")
    st.put_bytes(f"{rp}/uploads/sample_sheet.csv", sample_bytes, content_type="text/csv")

    uploaded = []
    for f in rsem_files:
        content = await f.read()
        name = Path(f.filename).name
        if not name:
            continue
        st.put_bytes(f"{rp}/uploads/{name}", content)
        uploaded.append(name)

    if not uploaded:
        raise HTTPException(status_code=400, detail="No RSEM files were uploaded")

    groups = _parse_groups_from_bytes(sample_bytes, {"sample_id", "file_name", "group"})

    meta.input_type = "rsem_expected_counts"
    meta.sample_sheet_name = sample_sheet.filename or "sample_sheet.csv"
    meta.count_matrix_name = None
    meta.uploaded_files = sorted(uploaded)
    meta.available_groups = groups
    meta.status = "uploaded"
    meta.error = None
    save_metadata(meta)

    return {
        "run_id": run_id,
        "status": meta.status,
        "input_type": meta.input_type,
        "uploaded_rsem_files": len(uploaded),
        "available_groups": groups,
    }


@app.post("/runs/{run_id}/upload-count-matrix")
async def upload_count_matrix_inputs(
    run_id: str,
    sample_sheet: UploadFile = File(...),
    count_matrix: UploadFile = File(...),
    input_type: InputType = Form(...),
) -> dict:
    try:
        meta = load_metadata(run_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    if input_type == "tpm_fpkm_only":
        raise HTTPException(
            status_code=400,
            detail="TPM/FPKM-only input is not valid for DESeq2. Please upload count-like data.",
        )

    if input_type not in COUNT_BASED_INPUT_TYPES - {"rsem_expected_counts"}:
        raise HTTPException(status_code=400, detail=f"Unsupported input_type for count matrix upload: {input_type}")

    st = get_storage()
    rp = run_prefix(run_id)

    sample_bytes = await sample_sheet.read()
    matrix_bytes = await count_matrix.read()
    if not sample_bytes or not matrix_bytes:
        raise HTTPException(status_code=400, detail="Sample sheet or count matrix is empty")

    st.put_bytes(f"{rp}/uploads/sample_sheet.csv", sample_bytes, content_type="text/csv")
    st.put_bytes(f"{rp}/uploads/count_matrix.csv", matrix_bytes, content_type="text/csv")

    groups = _parse_groups_from_bytes(sample_bytes, {"sample_id", "group"})

    meta.input_type = input_type
    meta.sample_sheet_name = sample_sheet.filename or "sample_sheet.csv"
    meta.count_matrix_name = count_matrix.filename or "count_matrix.csv"
    meta.uploaded_files = ["count_matrix.csv"]
    meta.available_groups = groups
    meta.status = "uploaded"
    meta.error = None
    save_metadata(meta)

    return {
        "run_id": run_id,
        "status": meta.status,
        "input_type": meta.input_type,
        "uploaded_count_matrix": True,
        "available_groups": groups,
    }


@app.post("/runs/{run_id}/start", response_model=StartRunResponse)
def start_run(run_id: str, req: StartRunRequest) -> StartRunResponse:
    try:
        meta = load_metadata(run_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    if meta.status not in {"uploaded", "failed", "completed"}:
        raise HTTPException(status_code=400, detail=f"Run status {meta.status} cannot be started")
    if req.max_labels < 0 or req.max_labels > 200:
        raise HTTPException(status_code=400, detail="max_labels must be between 0 and 200")

    if RUN_INLINE:
        process_run_job(
            run_id=run_id,
            group_a=req.group_a,
            group_b=req.group_b,
            label_significant_genes=req.label_significant_genes,
            max_labels=req.max_labels,
            run_gsea=req.run_gsea,
            prepare_cibersortx=req.prepare_cibersortx,
        )
        meta = load_metadata(run_id)
        return StartRunResponse(run_id=run_id, status=meta.status, job_id=None)

    q = get_queue()
    job = q.enqueue(
        "worker.jobs.process_run_job",
        run_id=run_id,
        group_a=req.group_a,
        group_b=req.group_b,
        label_significant_genes=req.label_significant_genes,
        max_labels=req.max_labels,
        run_gsea=req.run_gsea,
        prepare_cibersortx=req.prepare_cibersortx,
    )
    meta.status = "queued"
    meta.job_id = job.id
    meta.error = None
    save_metadata(meta)
    return StartRunResponse(run_id=run_id, status=meta.status, job_id=job.id)


@app.get("/runs/{run_id}")
def get_run(run_id: str) -> dict:
    try:
        meta = load_metadata(run_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    return meta.model_dump(mode="json")


@app.get("/runs/{run_id}/artifacts", response_model=ArtifactListResponse)
def get_artifacts(run_id: str) -> ArtifactListResponse:
    try:
        names = list_artifacts(run_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    st = get_storage()
    items = [ArtifactItem(name=n, size_bytes=st.file_size(artifact_key(run_id, n))) for n in names]
    return ArtifactListResponse(run_id=run_id, artifacts=items)


@app.get("/runs/{run_id}/artifacts/{artifact_name}")
def download_artifact(run_id: str, artifact_name: str):
    st = get_storage()
    key = artifact_key(run_id, artifact_name)
    if not st.exists(key):
        raise HTTPException(status_code=404, detail="Artifact not found")

    if st.mode == "s3":
        return RedirectResponse(url=st.presigned_url(key), status_code=307)

    path = run_dir(run_id) / "artifacts" / artifact_name
    return FileResponse(path)
