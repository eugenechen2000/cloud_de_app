from __future__ import annotations

from datetime import datetime

from pydantic import BaseModel


class CreateRunResponse(BaseModel):
    run_id: str
    project_id: str
    status: str
    created_at: datetime


class StartRunRequest(BaseModel):
    group_a: str | None = None
    group_b: str | None = None


class StartRunResponse(BaseModel):
    run_id: str
    status: str
    job_id: str | None


class ArtifactItem(BaseModel):
    name: str
    size_bytes: int


class ArtifactListResponse(BaseModel):
    run_id: str
    artifacts: list[ArtifactItem]
