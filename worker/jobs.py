from __future__ import annotations

from worker.pipeline import run_pipeline


def process_run_job(run_id: str, group_a: str | None = None, group_b: str | None = None) -> None:
    run_pipeline(run_id=run_id, group_a=group_a, group_b=group_b)
