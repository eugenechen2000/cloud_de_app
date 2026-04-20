from __future__ import annotations

from typing import Optional

from worker.pipeline import run_pipeline


def process_run_job(
    run_id: str,
    group_a: Optional[str] = None,
    group_b: Optional[str] = None,
    label_significant_genes: bool = False,
    max_labels: int = 25,
    run_gsea: bool = True,
    prepare_cibersortx: bool = True,
) -> None:
    run_pipeline(
        run_id=run_id,
        group_a=group_a,
        group_b=group_b,
        label_significant_genes=label_significant_genes,
        max_labels=max_labels,
        run_gsea=run_gsea,
        prepare_cibersortx=prepare_cibersortx,
    )
