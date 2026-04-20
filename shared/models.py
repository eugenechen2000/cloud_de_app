from __future__ import annotations

from datetime import datetime
from typing import Literal, Optional

from pydantic import BaseModel, Field

RunStatus = Literal["created", "uploaded", "queued", "running", "completed", "failed"]
InputType = Literal[
    "rsem_expected_counts",
    "raw_counts_matrix",
    "salmon_kallisto_estimated_counts",
    "tpm_fpkm_only",
]


class RunMetadata(BaseModel):
    run_id: str
    project_id: str
    status: RunStatus = "created"
    input_type: InputType = "rsem_expected_counts"
    created_at: datetime
    updated_at: datetime
    sample_sheet_name: Optional[str] = None
    count_matrix_name: Optional[str] = None
    uploaded_files: list[str] = Field(default_factory=list)
    available_groups: list[str] = Field(default_factory=list)
    selected_groups: list[str] = Field(default_factory=list)
    logs: list[str] = Field(default_factory=list)
    artifacts: list[str] = Field(default_factory=list)
    error: Optional[str] = None
    job_id: Optional[str] = None
