from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

from shared.config import DATA_ROOT
from shared.models import RunMetadata
from shared.storage import get_storage


def utcnow() -> datetime:
    return datetime.now(tz=timezone.utc)


def run_prefix(run_id: str) -> str:
    return f"runs/{run_id}"


def run_dir(run_id: str) -> Path:
    # Local working directory for current process.
    return DATA_ROOT / "runs" / run_id


def ensure_run_dirs(run_id: str) -> dict[str, Path]:
    base = run_dir(run_id)
    uploads = base / "uploads"
    work = base / "work"
    artifacts = base / "artifacts"
    for path in (base, uploads, work, artifacts):
        path.mkdir(parents=True, exist_ok=True)
    return {"base": base, "uploads": uploads, "work": work, "artifacts": artifacts}


def metadata_key(run_id: str) -> str:
    return f"{run_prefix(run_id)}/metadata.json"


def save_metadata(meta: RunMetadata) -> None:
    meta.updated_at = utcnow()
    payload = meta.model_dump_json(indent=2).encode("utf-8")
    get_storage().put_bytes(metadata_key(meta.run_id), payload, content_type="application/json")


def load_metadata(run_id: str) -> RunMetadata:
    st = get_storage()
    key = metadata_key(run_id)
    if not st.exists(key):
        raise FileNotFoundError(f"Run {run_id} not found")
    return RunMetadata.model_validate_json(st.get_bytes(key).decode("utf-8"))


def append_log(run_id: str, line: str) -> None:
    meta = load_metadata(run_id)
    meta.logs.append(line)
    save_metadata(meta)


def list_artifacts(run_id: str) -> list[str]:
    return get_storage().list_files(f"{run_prefix(run_id)}/artifacts")


def artifact_key(run_id: str, artifact_name: str) -> str:
    return f"{run_prefix(run_id)}/artifacts/{artifact_name}"
