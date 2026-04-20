from __future__ import annotations

from datetime import datetime
from typing import Optional

from pydantic import BaseModel


class CreateRunResponse(BaseModel):
    run_id: str
    project_id: str
    status: str
    created_at: datetime


class StartRunRequest(BaseModel):
    group_a: Optional[str] = None
    group_b: Optional[str] = None
    label_significant_genes: bool = False
    max_labels: int = 25
    run_gsea: bool = True
    prepare_cibersortx: bool = True


class StartRunResponse(BaseModel):
    run_id: str
    status: str
    job_id: Optional[str]


class ArtifactItem(BaseModel):
    name: str
    size_bytes: int


class ArtifactListResponse(BaseModel):
    run_id: str
    artifacts: list[ArtifactItem]
